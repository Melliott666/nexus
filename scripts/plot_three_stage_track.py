#!/usr/bin/env python3
"""Compare a KingCRAB truth trail, EL charge image, and II optical image."""

from argparse import ArgumentParser, SUPPRESS
from pathlib import Path
import csv
import shutil
import subprocess
import sys
import tempfile


def decode(values):
    import numpy as np

    if values.dtype.kind in ("S", "a"):
        return np.char.decode(values, "utf-8", errors="ignore")
    return values.astype(str)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("h5_file", type=Path)
    parser.add_argument("--event", type=int, default=0)
    parser.add_argument(
        "--all-events", action="store_true",
        help="analyze every event and bundle the results into one ZIP archive",
    )
    parser.add_argument(
        "--archive", type=Path, default=None,
        help="output ZIP for --all-events (default: <input>_three_stage_tracks.zip)",
    )
    parser.add_argument("--z-plane", type=float, default=1449.4)
    parser.add_argument("--ii-center-x", type=float, default=-82.55)
    parser.add_argument("--ii-center-y", type=float, default=143.002)
    parser.add_argument("--ii-radius", type=float, default=25.4)
    parser.add_argument("--slice-z-min", type=float, default=1430.0)
    parser.add_argument("--slice-z-max", type=float, default=1450.0)
    parser.add_argument("--slice-z-step", type=float, default=0.25)
    parser.add_argument("--diagnostic-radius", type=float, default=100.0)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--skip-stack-zip", action="store_true", help=SUPPRESS)
    args = parser.parse_args()
    if args.slice_z_step <= 0 or args.slice_z_max < args.slice_z_min:
        parser.error("Require slice-z-step > 0 and slice-z-max >= slice-z-min")

    import h5py
    import matplotlib
    import numpy as np

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    h5_file = args.h5_file.expanduser().resolve()

    if args.all_events:
        with h5py.File(h5_file, "r") as h5:
            for table in ("/MC/hits", "/MC/particles", "/DEBUG/steps"):
                if table not in h5:
                    raise KeyError(f"Missing {table}")
            event_ids = np.unique(h5["/MC/hits"]["event_id"]).astype(int)
        if len(event_ids) == 0:
            raise RuntimeError("No events found in /MC/hits")

        archive = args.archive
        if archive is None:
            archive = h5_file.with_name(
                f"{h5_file.stem}_three_stage_tracks.zip"
            )
        archive = archive.expanduser().resolve()
        if archive.suffix.lower() != ".zip":
            archive = archive.with_suffix(".zip")
        archive.parent.mkdir(parents=True, exist_ok=True)

        with tempfile.TemporaryDirectory(prefix="kingcrab_muon_events_") as tmp:
            bundle = Path(tmp) / f"{h5_file.stem}_three_stage_tracks"
            bundle.mkdir()
            for event_id in event_ids:
                event_dir = bundle / f"event_{event_id:06d}"
                event_dir.mkdir()
                event_output = event_dir / "three_stage_track.png"
                command = [
                    sys.executable, str(Path(__file__).resolve()), str(h5_file),
                    "--event", str(event_id),
                    "--z-plane", str(args.z_plane),
                    "--ii-center-x", str(args.ii_center_x),
                    "--ii-center-y", str(args.ii_center_y),
                    "--ii-radius", str(args.ii_radius),
                    "--slice-z-min", str(args.slice_z_min),
                    "--slice-z-max", str(args.slice_z_max),
                    "--slice-z-step", str(args.slice_z_step),
                    "--diagnostic-radius", str(args.diagnostic_radius),
                    "--output", str(event_output),
                    "--skip-stack-zip",
                ]
                print(f"Analyzing event {event_id} ({event_dir.name})", flush=True)
                subprocess.run(command, check=True)

            manifest = bundle / "README.txt"
            manifest.write_text(
                f"Source HDF5 file: {h5_file}\n"
                f"Events: {', '.join(map(str, event_ids))}\n"
                f"Number of events: {len(event_ids)}\n"
                "Each event_NNNNNN directory contains that event's plots, "
                "summary, and focal-plane image stack.\n"
            )
            archive_base = archive.with_suffix("")
            created = Path(shutil.make_archive(
                str(archive_base), "zip", root_dir=bundle.parent,
                base_dir=bundle.name,
            ))
        print(f"Bundled {len(event_ids)} events into {created}")
        return

    output = args.output
    if output is None:
        output = h5_file.with_name(f"{h5_file.stem}_three_stage_track.png")
    output = output.expanduser().resolve()
    summary_path = output.with_suffix(".txt")

    with h5py.File(h5_file, "r") as h5:
        for table in ("/MC/hits", "/MC/particles", "/DEBUG/steps"):
            if table not in h5:
                raise KeyError(f"Missing {table}")
        hits = h5["/MC/hits"][:]
        particles = h5["/MC/particles"][:]
        steps = h5["/DEBUG/steps"][:]

    hits = hits[hits["event_id"] == args.event]
    particles = particles[particles["event_id"] == args.event]
    steps = steps[steps["event_id"] == args.event]
    if len(hits) == 0:
        raise RuntimeError(f"No /MC/hits rows for event {args.event}")

    step_particle = decode(steps["particle_name"])
    initial_volume = decode(steps["initial_volume"])
    final_volume = decode(steps["final_volume"])
    process = decode(steps["proc_name"])

    entry_mask = np.char.startswith(process, "EL_ELECTRON_ENTRY:")
    entry_steps = steps[entry_mask]
    if len(entry_steps) == 0:
        raise RuntimeError(
            "No EL_ELECTRON_ENTRY rows. Run with record_el_electron_entries true."
        )
    # The fast drift's post-step x-y is the diffused charge-arrival position.
    entry_x = entry_steps["final_x"].astype(float)
    entry_y = entry_steps["final_y"].astype(float)

    scan_mask = (
        (step_particle == "opticalphoton")
        & (initial_volume == "FOCAL_SCAN")
        & (steps["initial_z"] != steps["final_z"])
    )
    scan = steps[scan_mask]
    def crossings_at(z_plane):
        z0 = scan["initial_z"].astype(float)
        z1 = scan["final_z"].astype(float)
        crossing = ((z0 <= z_plane) & (z1 >= z_plane)) | (
            (z1 <= z_plane) & (z0 >= z_plane)
        )
        selected = scan[crossing]
        if len(selected) == 0:
            return np.array([]), np.array([])
        fraction = (z_plane - selected["initial_z"]) / (
            selected["final_z"] - selected["initial_z"]
        )
        x = selected["initial_x"] + fraction * (
            selected["final_x"] - selected["initial_x"]
        )
        y = selected["initial_y"] + fraction * (
            selected["final_y"] - selected["initial_y"]
        )
        image_keys = np.rec.fromarrays(
            [selected["event_id"], selected["particle_id"]],
            names=("event_id", "particle_id"),
        )
        _, unique = np.unique(image_keys, return_index=True)
        return np.asarray(x)[unique], np.asarray(y)[unique]

    image_x, image_y = crossings_at(args.z_plane)
    if len(image_x) == 0:
        raise RuntimeError(f"No photons cross z={args.z_plane:g} mm")
    inside_ii = np.hypot(
        image_x - args.ii_center_x, image_y - args.ii_center_y
    ) <= args.ii_radius
    image_x = image_x[inside_ii] - args.ii_center_x
    image_y = image_y[inside_ii] - args.ii_center_y

    primary = particles[particles["primary"] != 0]
    particle_names = decode(primary["particle_name"]) if len(primary) else np.array([])
    muon = primary[np.char.find(particle_names, "mu") >= 0] if len(primary) else primary

    # Preserve the full three-dimensional truth information in one convenient
    # file. Each panel contains the same /MC/hits data, viewed from a different
    # camera angle; no trajectory is inferred from optical-photon endpoints.
    truth_3d_output = output.with_name(f"{output.stem}_truth_3d_views.png")
    fig_3d = plt.figure(figsize=(19, 6), constrained_layout=True)
    camera_views = (
        (22, -60, "Oblique view"),
        (22, 30, "Opposite oblique view"),
        (75, -90, "Near top view"),
    )
    truth_3d = None
    hit_x = hits["x"].astype(float)
    hit_y = hits["y"].astype(float)
    hit_z = hits["z"].astype(float)
    spans = np.ptp(np.column_stack((hit_x, hit_y, hit_z)), axis=0)
    spans = np.maximum(spans, 1.0)
    for panel, (elevation, azimuth, title) in enumerate(camera_views, start=1):
        axis = fig_3d.add_subplot(1, 3, panel, projection="3d")
        truth_3d = axis.scatter(
            hit_x, hit_y, hit_z, c=hits["energy"], s=10,
            cmap="inferno", linewidths=0, depthshade=False,
        )
        if len(muon):
            row = muon[0]
            axis.plot(
                [row["initial_x"], row["final_x"]],
                [row["initial_y"], row["final_y"]],
                [row["initial_z"], row["final_z"]],
                color="tab:cyan", linewidth=1.4,
                label="primary muon endpoints",
            )
            axis.legend(loc="upper right", fontsize=7)
        axis.view_init(elev=elevation, azim=azimuth)
        axis.set_box_aspect(spans)
        axis.set(
            xlabel="x [mm]", ylabel="y [mm]", zlabel="z [mm]",
            title=f"{title}\nelev={elevation}°, azim={azimuth}°",
        )
        axis.grid(alpha=0.25)
    fig_3d.colorbar(
        truth_3d, ax=fig_3d.axes, shrink=0.72, pad=0.04,
        label="energy deposit [MeV]",
    )
    fig_3d.suptitle(
        f"KingCRAB event {args.event}: original muon /MC/hits trail",
        fontsize=15,
    )
    truth_3d_output.parent.mkdir(parents=True, exist_ok=True)
    fig_3d.savefig(truth_3d_output, dpi=210, bbox_inches="tight")
    plt.close(fig_3d)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)

    truth = axes[0].scatter(
        hits["x"], hits["y"], c=hits["energy"], s=15,
        cmap="inferno", linewidths=0,
    )
    if len(muon):
        row = muon[0]
        axes[0].plot(
            [row["initial_x"], row["final_x"]],
            [row["initial_y"], row["final_y"]],
            color="tab:cyan", linewidth=1.3, label="primary muon endpoints",
        )
        axes[0].legend(loc="best", fontsize=8)
    fig.colorbar(truth, ax=axes[0], label="energy deposit [MeV]")
    axes[0].set_title(f"1. Truth energy-deposition trail\n{len(hits)} hits")

    el_hist = axes[1].hist2d(entry_x, entry_y, bins=100, cmap="viridis")
    fig.colorbar(el_hist[3], ax=axes[1], label="arriving ionization electrons/bin")
    axes[1].set_title(f"2. Charge arriving at EL gap\n{len(entry_x)} electrons")

    axes[2].scatter(image_x, image_y, s=4, alpha=0.45, linewidths=0)
    axes[2].add_patch(plt.Circle(
        (0, 0), args.ii_radius, fill=False, color="tab:red", linewidth=1.4,
        label="physical II aperture",
    ))
    axes[2].legend(loc="best", fontsize=8)
    axes[2].set(
        xlim=(-args.ii_radius, args.ii_radius),
        ylim=(-args.ii_radius, args.ii_radius),
        title=(
            f"3. Optical image at z={args.z_plane:.1f} mm\n"
            f"{len(image_x)} photons inside II aperture"
        ),
        xlabel="II-local x [mm]", ylabel="II-local y [mm]",
    )

    for axis in axes[:2]:
        axis.set(xlabel="x [mm]", ylabel="y [mm]")
    for axis in axes:
        axis.set_aspect("equal", adjustable="box")
        axis.grid(alpha=0.2)
    fig.suptitle(
        f"KingCRAB event {args.event}: truth → drifted charge → optical image",
        fontsize=15,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=210, bbox_inches="tight")
    plt.close(fig)

    deposited_mev = float(np.sum(hits["energy"]))
    summary = (
        f"HDF5 file: {h5_file}\n"
        f"Event: {args.event}\n"
        f"Truth energy-deposition hits: {len(hits)}\n"
        f"Total deposited energy: {deposited_mev:.6g} MeV\n"
        f"Ionization electrons recorded on first EL entry: {len(entry_x)}\n"
        f"Photons crossing z={args.z_plane:.6f} mm inside II aperture: {len(image_x)}\n"
    )
    summary_path.write_text(summary)
    print(summary, end="")
    print(f"Saved {output}")
    print(f"Saved {truth_3d_output}")
    print(f"Saved {summary_path}")

    # Save a conventional stack of independent 2D optical images so the
    # track's focus and distortion can be inspected plane by plane.
    stack_dir = output.with_name(f"{output.stem}_ii_focal_stack")
    stack_dir.mkdir(parents=True, exist_ok=True)
    z_planes = np.arange(
        args.slice_z_min,
        args.slice_z_max + 0.5 * args.slice_z_step,
        args.slice_z_step,
    )
    stack_counts = []
    for z_plane in z_planes:
        plane_x, plane_y = crossings_at(z_plane)
        local_x = plane_x - args.ii_center_x
        local_y = plane_y - args.ii_center_y
        aperture = np.hypot(local_x, local_y) <= args.ii_radius
        stack_counts.append((z_plane, len(local_x), int(np.count_nonzero(aperture))))

        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(
            local_x[~aperture], local_y[~aperture], s=3, alpha=0.25,
            color="0.55", linewidths=0, label="outside physical II aperture",
        )
        ax.scatter(
            local_x[aperture], local_y[aperture], s=4, alpha=0.5,
            color="tab:blue", linewidths=0, label="inside physical II aperture",
        )
        ax.add_patch(plt.Circle(
            (0, 0), args.ii_radius, fill=False, color="tab:red", linewidth=1.5,
        ))
        ax.set(
            xlim=(-args.diagnostic_radius, args.diagnostic_radius),
            ylim=(-args.diagnostic_radius, args.diagnostic_radius),
            xlabel="II-local x [mm]", ylabel="II-local y [mm]",
            title=(
                f"Muon optical track at z = {z_plane:.3f} mm\n"
                f"wide cylinder = {len(local_x)}; inside II = {np.count_nonzero(aperture)}"
            ),
        )
        ax.legend(loc="upper right", fontsize=8)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(alpha=0.2)
        plane_label = f"{z_plane:08.3f}".replace("-", "m").replace(".", "p")
        filename = f"ii_track_z_{plane_label}_mm.png"
        fig.savefig(stack_dir / filename, dpi=180, bbox_inches="tight")
        plt.close(fig)

    with (stack_dir / "ii_focal_stack_counts.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["z_mm", "wide_cylinder_photons", "inside_ii_photons"])
        writer.writerows(stack_counts)
    print(f"Saved {len(z_planes)} II focal-plane images to {stack_dir}")
    if not args.skip_stack_zip:
        archive = shutil.make_archive(str(stack_dir), "zip", stack_dir)
        print(f"Created {archive}")


if __name__ == "__main__":
    main()
