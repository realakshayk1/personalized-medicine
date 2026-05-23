import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "Lattice",
  description: "AI-powered single-cell RNA-seq analysis",
};

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html lang="en" className="dark h-full">
      <head>
        <meta name="color-scheme" content="dark" />
      </head>
      <body className="h-full overflow-hidden">{children}</body>
    </html>
  );
}
