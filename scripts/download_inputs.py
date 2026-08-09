from pathlib import Path
import argparse
from urllib.request import urlopen, Request

def download(url, out):
    out=Path(out); out.parent.mkdir(parents=True,exist_ok=True)
    req=Request(url,headers={'User-Agent':'ICH-expression-contrast-reproducibility/1.0'})
    with urlopen(req,timeout=120) as src, out.open('wb') as dst:
        while True:
            block=src.read(1024*1024)
            if not block: break
            dst.write(block)

p=argparse.ArgumentParser()
p.add_argument('--expr-url',required=True); p.add_argument('--expr-out',required=True)
p.add_argument('--gmt-url',required=True); p.add_argument('--gmt-out',required=True)
a=p.parse_args()
download(a.expr_url,a.expr_out)
download(a.gmt_url,a.gmt_out)

