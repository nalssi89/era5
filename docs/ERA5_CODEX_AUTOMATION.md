# ERA5 Codex Automation

This repository records the operating contract for the daily ERA5 dashboard
automation. The schedule itself is managed by Codex automation under the
`era5` automation id; this file is the versioned reference for the run.

## Daily Run

The automation runs once per day at 03:00 UTC, which is 12:00 in Korea
(Asia/Seoul). It uses the local workspace `D:\WORK\Projects\Codex\ERA5`.

1. Read `outputs/wave_dashboard/summary.json` and the current public dashboard.
2. Run the latest-data update with a five-day ERA5T lag:

   ```powershell
   powershell.exe -NoProfile -ExecutionPolicy Bypass -File ops\update_era5_dashboard_daily.ps1 -Start 2026-04-01 -End auto -LagDays 5
   ```

3. Run the strict local quality gate:

   ```powershell
   py -3 scripts\qa_wave_dashboard.py --strict
   ```

4. Compare the local latest date and dashboard files with the public site.
5. Synchronize and publish only when the latest date or dashboard artifacts
   changed. Do not redeploy an unchanged dashboard.

## Publication Contract

The public dashboard is published at:

`https://era5-korea-climate-monitor.solverrrrr.chatgpt.site/dashboard/index.html`

The deployment must preserve the dashboard's daily, Monday-to-Sunday weekly,
and calendar-month views. The public response and the latest PNG are checked
after deployment. Unrelated `drizzle`, `examples`, and `tests` changes are not
part of the dashboard publication commit.

## Repository Boundaries

The Codex scheduler and Sites deployment credentials are external to GitHub.
They must not be copied into this repository. CDS credentials, source-repository
tokens, and deployment tokens remain in their respective secret stores.

The machine-local implementation lives in the ERA5 workspace referenced above;
this repository keeps the public dashboard assets and the versioned automation
contract so the process is auditable and reproducible.
