# Packs PNGs into an .ico (PNG entries, accepted since Windows Vista). Standard library only.
# Uso: python3 -I make_ico.py saida.ico 16.png 32.png ... (cada PNG quadrado, ≤ 256 px)
import struct
import sys

out, pngs = sys.argv[1], sys.argv[2:]
blobs = [open(p, "rb").read() for p in pngs]
head = struct.pack("<HHH", 0, 1, len(blobs))
off = 6 + 16 * len(blobs)
dirs = b""
for b in blobs:
    w, h = struct.unpack(">II", b[16:24])  # IHDR
    dirs += struct.pack("<BBBBHHII", w % 256, h % 256, 0, 0, 1, 32, len(b), off)
    off += len(b)
open(out, "wb").write(head + dirs + b"".join(blobs))
