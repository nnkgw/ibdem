#!/usr/bin/env python3
"""Generate animation.gif from capture_fr*.bmp files.

Usage:
  python make_gif.py [<bmp_dir>] [--out <path>] [--fps <n>] [--scale <factor>]

Defaults: bmp_dir="./"  out="animation.gif"  fps=5  scale=0.5
"""
import glob, sys, os, argparse
from PIL import Image

parser = argparse.ArgumentParser()
parser.add_argument("bmp_dir", nargs="?", default=".")
parser.add_argument("--out",   default="animation.gif")
parser.add_argument("--fps",   type=float, default=5.0)
parser.add_argument("--scale", type=float, default=0.5)
args = parser.parse_args()

pattern = os.path.join(args.bmp_dir, "capture_fr*.bmp")
files   = sorted(glob.glob(pattern))
if not files:
    print(f"No files matching {pattern}", file=sys.stderr)
    sys.exit(1)

duration_ms = int(1000 / args.fps)

frames = []
for f in files:
    img = Image.open(f).convert("RGB")
    if args.scale != 1.0:
        w, h = img.size
        img = img.resize((int(w * args.scale), int(h * args.scale)), Image.LANCZOS)
    frames.append(img.convert("P", palette=Image.ADAPTIVE, colors=128))

frames[0].save(
    args.out,
    save_all=True,
    append_images=frames[1:],
    duration=duration_ms,
    loop=0,
    optimize=True,
)
print(f"Saved {args.out}  ({len(frames)} frames, {os.path.getsize(args.out)//1024} KB)")
