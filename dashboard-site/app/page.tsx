import type { Metadata } from "next";

export const metadata: Metadata = {
  title: "ERA5 Korea Climate Monitor",
  description: "ERA5 재분석장 기후감시 대시보드",
};

export default function Home() {
  return (
    <main className="dashboard-shell">
      <iframe
        className="dashboard-frame"
        src="/dashboard/index.html"
        title="ERA5 Korea Climate Monitor"
      />
    </main>
  );
}
