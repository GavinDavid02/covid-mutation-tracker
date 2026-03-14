import zstandard as zstd, io

with open('data/raw/metadata.tsv.zst', 'rb') as f:
    dctx = zstd.ZstdDecompressor()
    stream = dctx.stream_reader(f)
    reader = io.TextIOWrapper(stream, encoding='utf-8')
    header = reader.readline().strip().split('\t')
    row1 = reader.readline().strip().split('\t')

print('COLUMNS:')
for i, (h, v) in enumerate(zip(header, row1)):
    print(f'  {i:3d}  {h:40s}  {v[:60]}')
