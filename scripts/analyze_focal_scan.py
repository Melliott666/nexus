#!/usr/bin/env python3
"""Reconstruct KingCRAB optical images on virtual planes in FOCAL_SCAN."""

from argparse import ArgumentParser
from pathlib import Path
import csv
import shutil


def decode(values):
    import numpy as np

    if values.dtype.kind in ("S", "a"):
        return np.char.decode(values, "utf-8", errors="ignore")
    return values.astype(str)


def plane_name(z):
    return f"focal_scan_z_{z:08.3f}_mm".replace("-", "m").replace(".", "p")


def containment_radius(x, y, percentile):
    import numpy as np

    if len(x) == 0:
        return float("nan")
    radius = np.hypot(x - np.mean(x), y - np.mean(y))
    return float(np.percentile(radius, percentile))


def main():
    parser = ArgumentParser(
        description=(
            "Create virtual image-plane plots from optical-photon segments "
            "stored by SaveAllSteppingAction in /DEBUG/steps."
        )
    )
    parser.add_argument("h5_file", type=Path)
    parser.add_argument("--z-min", type=float, default=1447.0)
    parser.add_argument("--z-max", type=float, default=1451.0)
    parser.add_argument("--z-step", type=float, default=0.05)
    parser.add_argument("--center-x", type=float, default=-82.55)
    parser.add_argument("--center-y", type=float, default=143.002)
    parser.add_argument("--radius", type=float, default=30.0)
    parser.add_argument(
        "--roi-radius", type=float, default=4.0,
        help=(
            "Radius in mm around the designed optical axis used to isolate "
            "the compact, desired alpha image. Default: 4.0"
        ),
    )
    parser.add_argument("--volume", default="FOCAL_SCAN")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument(
        "--no-zip", action="store_true", help="Do not create a ZIP archive."
    )
    args = parser.parse_args()

    if args.z_step <= 0:
        parser.error("--z-step must be positive")
    if args.z_max < args.z_min:
        parser.error("--z-max must be greater than or equal to --z-min")

    import h5py
    import matplotlib
    import numpy as np

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    h5_file = args.h5_file.expanduser().resolve()
    if not h5_file.exists():
        raise FileNotFoundError(h5_file)

    output_dir = args.output_dir
    if output_dir is None:
        output_dir = h5_file.with_name(f"{h5_file.stem}_focal_scan")
    output_dir = output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    with h5py.File(h5_file, "r") as h5:
        if "/DEBUG/steps" not in h5:
            raise KeyError(
                "Missing /DEBUG/steps. Run with SaveAllSteppingAction and "
                "select the FOCAL_SCAN volume."
            )

        steps = h5["/DEBUG/steps"][:]

    required = {
        "event_id", "particle_id", "particle_name", "initial_volume",
        "final_volume", "initial_x", "initial_y", "initial_z", "final_x",
        "final_y", "final_z",
    }
    missing = required.difference(steps.dtype.names or ())
    if missing:
        raise KeyError(f"/DEBUG/steps is missing: {', '.join(sorted(missing))}")

    particle = decode(steps["particle_name"])
    initial_volume = decode(steps["initial_volume"])
    final_volume = decode(steps["final_volume"])

    # Retain segments whose starting point is inside the scan volume. Entry
    # steps end at the boundary and contain no path through the cylinder;
    # internal and exit steps begin inside and span its usable ray trajectory.
    mask = (
        (particle == "opticalphoton")
        & (np.char.find(initial_volume, args.volume) >= 0)
        & (steps["final_z"] != steps["initial_z"])
    )
    segments = steps[mask]
    if len(segments) == 0:
        raise RuntimeError(f"No optical-photon segments begin in {args.volume}")

    z_planes = np.arange(
        args.z_min, args.z_max + 0.5 * args.z_step, args.z_step
    )
    metrics = []

    for z_plane in z_planes:
        z0 = segments["initial_z"].astype(float)
        z1 = segments["final_z"].astype(float)
        crosses = ((z0 <= z_plane) & (z1 >= z_plane)) | (
            (z0 >= z_plane) & (z1 <= z_plane)
        )
        selected = segments[crosses]

        if len(selected):
            z0 = selected["initial_z"].astype(float)
            z1 = selected["final_z"].astype(float)
            fraction = (z_plane - z0) / (z1 - z0)
            x = selected["initial_x"] + fraction * (
                selected["final_x"] - selected["initial_x"]
            )
            y = selected["initial_y"] + fraction * (
                selected["final_y"] - selected["initial_y"]
            )

            # A plane exactly on a step boundary can select two adjacent steps.
            # Keep only one crossing for each event/track pair.
            keys = np.rec.fromarrays(
                [selected["event_id"], selected["particle_id"]],
                names=("event_id", "particle_id"),
            )
            _, unique_indices = np.unique(keys, return_index=True)
            x = np.asarray(x)[unique_indices]
            y = np.asarray(y)[unique_indices]
        else:
            x = np.array([])
            y = np.array([])

        mean_x = float(np.mean(x)) if len(x) else float("nan")
        mean_y = float(np.mean(y)) if len(y) else float("nan")
        rms_radius = (
            float(np.sqrt(np.mean((x - mean_x) ** 2 + (y - mean_y) ** 2)))
            if len(x) else float("nan")
        )
        r68 = containment_radius(x, y, 68.0)
        r95 = containment_radius(x, y, 95.0)

        # Analyze the compact image independently of the displaced curved
        # ghost. The ROI is fixed around the designed image-axis location so a
        # changing ghost population cannot pull the selection centroid.
        distance_from_axis = np.hypot(
            x - args.center_x, y - args.center_y
        )
        roi = distance_from_axis < args.roi_radius
        roi_x = x[roi]
        roi_y = y[roi]
        roi_mean_x = float(np.mean(roi_x)) if len(roi_x) else float("nan")
        roi_mean_y = float(np.mean(roi_y)) if len(roi_y) else float("nan")
        roi_rms = (
            float(np.sqrt(np.mean(
                (roi_x - roi_mean_x) ** 2 + (roi_y - roi_mean_y) ** 2
            ))) if len(roi_x) else float("nan")
        )
        roi_r68 = containment_radius(roi_x, roi_y, 68.0)
        roi_r95 = containment_radius(roi_x, roi_y, 95.0)

        metrics.append((
            z_plane, len(x), mean_x, mean_y, rms_radius, r68, r95,
            len(roi_x), roi_mean_x, roi_mean_y, roi_rms, roi_r68, roi_r95,
        ))

        fig, ax = plt.subplots(figsize=(7, 7))
        ax.scatter(x, y, s=3, alpha=0.45, linewidths=0)
        ax.set(
            xlim=(args.center_x - args.radius, args.center_x + args.radius),
            ylim=(args.center_y - args.radius, args.center_y + args.radius),
            xlabel="x [mm]",
            ylabel="y [mm]",
            title=(
                f"KingCRAB virtual image at z = {z_plane:.3f} mm\n"
                f"all N = {len(x)}; compact N = {len(roi_x)}, "
                f"compact r68 = {roi_r68:.3f} mm"
            ),
        )
        roi_circle = plt.Circle(
            (args.center_x, args.center_y), args.roi_radius,
            fill=False, color="tab:red", linestyle="--", linewidth=1.2,
            label="compact-image ROI",
        )
        ax.add_patch(roi_circle)
        ax.legend(loc="lower right")
        ax.set_aspect("equal", adjustable="box")
        ax.grid(alpha=0.25)
        fig.savefig(output_dir / f"{plane_name(z_plane)}.png", dpi=180,
                    bbox_inches="tight")
        plt.close(fig)

    metrics_path = output_dir / "focal_scan_metrics.csv"
    with metrics_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["z_mm", "photons", "mean_x_mm", "mean_y_mm", "rms_radius_mm",
             "r68_mm", "r95_mm", "roi_photons", "roi_mean_x_mm",
             "roi_mean_y_mm", "roi_rms_radius_mm", "roi_r68_mm",
             "roi_r95_mm"]
        )
        writer.writerows(metrics)

    # Column 11 is the compact-image r68. This, rather than the all-photon
    # r68 in column 5, determines the desired image-intensifier location.
    finite_metrics = [row for row in metrics if np.isfinite(row[11])]
    if finite_metrics:
        best = min(finite_metrics, key=lambda row: row[11])
        summary = (
            f"Minimum compact-image r68 plane: z = {best[0]:.6f} mm\n"
            f"Compact-image photons: {best[7]}\n"
            f"Compact-image centroid: ({best[8]:.6f}, {best[9]:.6f}) mm\n"
            f"Compact-image RMS radius: {best[10]:.6f} mm\n"
            f"Compact-image r68: {best[11]:.6f} mm\n"
            f"Compact-image r95: {best[12]:.6f} mm\n"
            f"ROI radius about designed axis: {args.roi_radius:.6f} mm\n"
        )
    else:
        summary = "No photon crossings were found on the requested planes.\n"
    (output_dir / "focal_scan_summary.txt").write_text(summary)
    print(summary, end="")
    print(f"Saved {len(z_planes)} plane images and metrics to {output_dir}")

    if not args.no_zip:
        archive = shutil.make_archive(str(output_dir), "zip", output_dir)
        print(f"Created {archive}")


if __name__ == "__main__":
    main()
