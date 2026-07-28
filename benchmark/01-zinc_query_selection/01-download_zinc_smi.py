import asyncio
import sys
from pathlib import Path

URL_FILE = Path(__file__).parent / "tranche_urls.uri"
OUTPUT_DIR = Path(__file__).parent / "../data/zinc_tranches"
FAILED_FILE = Path(__file__).parent / "failed_urls.uri"
CONCURRENCY = 5  # max 5 concurrent downloads set by ZINC server


async def download(url: str, dest: Path) -> bool:
    dest.parent.mkdir(parents=True, exist_ok=True)
    proc = await asyncio.create_subprocess_exec(
        "curl", "--remote-time", "--fail", "--silent", "--show-error",
        "--create-dirs", "-o", str(dest), url,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    _, stderr = await proc.communicate()
    if proc.returncode == 0:
        print(f"Downloaded: {dest.name}")
        return True
    else:
        print(f"FAILED: {url} — {stderr.decode().strip()}")
        if dest.exists():
            dest.unlink()
        return False


async def main() -> None:
    urls = [line.strip() for line in URL_FILE.read_text().splitlines() if line.strip()]
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # skip already-downloaded files to support resuming
    pending = [u for u in urls if not (OUTPUT_DIR / u.split("/")[-1]).exists()]
    already_done = len(urls) - len(pending)
    if already_done:
        print(f"Skipping {already_done} already-downloaded file(s).")
    print(f"Downloading {len(pending)} file(s)...")

    semaphore = asyncio.Semaphore(CONCURRENCY)
    failed: list[str] = []
    lock = asyncio.Lock()

    async def bounded_download(url: str) -> None:
        dest = OUTPUT_DIR / url.split("/")[-1]
        async with semaphore:
            ok = await download(url, dest)
            if not ok:
                async with lock:
                    failed.append(url)

    await asyncio.gather(*[bounded_download(url) for url in pending])

    if failed:
        FAILED_FILE.write_text("\n".join(failed) + "\n")
        print(f"\n{len(failed)} file(s) failed. URLs saved to {FAILED_FILE.name}")
        print(f"To retry, run:  python {Path(__file__).name} --resume-failed")
    else:
        if FAILED_FILE.exists():
            FAILED_FILE.unlink()
        print(f"\nDone. All {len(pending)} file(s) downloaded to {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    if "--resume-failed" in sys.argv:
        if not FAILED_FILE.exists():
            print("No failed_urls.uri found — nothing to resume.")
            sys.exit(0)
        URL_FILE = FAILED_FILE  # noqa: F811
    asyncio.run(main())
