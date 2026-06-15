# ras-commander-hydro — Agent Notes

## Documentation Site

Docs publish to **https://rascommander.info/hydro** on every push to `main` (self-hosted; build infra
in `CLB-Engineering-Corporation/ras-commander-docs`). A broken `mkdocs.yml` fails the live build.

- This is a **tool** exposing a focused subset of the ras-commander library. The cohesion hub injects
  a "limited subset — for the full library, use ras-commander" banner from the shared theme; keep that
  framing and point users to **https://rascommander.info/ras** for the full library.
- Mechanics-forward authoring voice: defer methodology, parameter appropriateness, and
  regulatory/standard-of-care questions to HEC's manuals and the reader's regional/agency references.
