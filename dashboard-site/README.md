# ERA5 Korea Climate Monitor

This directory contains the externally published ERA5 dashboard source.

Live site: https://era5-korea-climate-monitor.solverrrrr.chatgpt.site

The dashboard includes the latest validated daily, weekly, and monthly ERA5
fields plus WAF and RWS proxy diagnostics. The `public/dashboard/` directory
contains the self-contained dashboard HTML and its PNG figure assets.

## Build

```bash
npm install
npm run build
```

The site shell serves `public/dashboard/index.html` at the `/dashboard/` path.
