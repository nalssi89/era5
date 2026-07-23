import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "ERA5 Korea Climate Monitor",
  description: "ERA5 재분석장 기후감시 대시보드",
  icons: {
    icon: "/favicon.svg",
    shortcut: "/favicon.svg",
  },
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="ko">
      <body>{children}</body>
    </html>
  );
}
