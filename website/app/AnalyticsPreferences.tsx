"use client";

import { useEffect, useState } from "react";

type Choice = "allow" | "decline";
const preference = "fixpoint:google-analytics";
const measurementId = "G-9279HK4C5Y";
const disableKey = "ga-disable-G-9279HK4C5Y";

declare global {
  interface Window {
    dataLayer?: IArguments[];
    [disableKey]?: boolean;
  }
}

function startAnalytics() {
  if (document.getElementById("fixpoint-google-tag")) return;
  window[disableKey] = false;
  const queue = window.dataLayer ??= [];
  function gtag(..._args: unknown[]) { queue.push(arguments); }
  gtag("consent", "default", {
    analytics_storage: "granted", ad_storage: "denied",
    ad_user_data: "denied", ad_personalization: "denied",
  });
  gtag("js", new Date());
  gtag("config", measurementId, {
    allow_google_signals: false,
    allow_ad_personalization_signals: false,
    cookie_domain: location.hostname,
  });
  const script = document.createElement("script");
  script.id = "fixpoint-google-tag";
  script.async = true;
  script.src = "https://www.googletagmanager.com/gtag/js?id=" + measurementId;
  document.head.appendChild(script);
}

export default function AnalyticsPreferences() {
  const [enabled, setEnabled] = useState(false);
  const [open, setOpen] = useState(false);
  const [choice, setChoice] = useState<Choice | null>(null);

  useEffect(() => {
    if (!["fixpoint-theory.com", "www.fixpoint-theory.com", "prime.fixpoint-theory.com"].includes(location.hostname)) return;
    let saved: string | null = null;
    try { saved = localStorage.getItem(preference); } catch { /* Storage is optional. */ }
    if (saved === "allow" || saved === "decline") setChoice(saved);
    else setOpen(true);
    if (saved === "allow") startAnalytics();
    setEnabled(true);
  }, []);

  function choose(next: Choice) {
    try { localStorage.setItem(preference, next); } catch { /* Apply the choice for this visit. */ }
    setOpen(false);
    setChoice(next);
    if (next === "allow") {
      startAnalytics();
    } else if (choice === "allow") {
      window[disableKey] = true;
      for (const pair of document.cookie.split(";")) {
        const name = pair.trim().split("=")[0];
        if (name === "_ga" || name.startsWith("_ga_")) {
          document.cookie = name + "=; Max-Age=0; Path=/; Domain=" + location.hostname;
        }
      }
      location.reload();
    }
  }

  if (!enabled) return null;
  return (
    <div className="px-4 pb-4 text-center text-sm text-slate-700">
      <button type="button" className="underline" onClick={() => setOpen(true)}>Analytics preferences</button>
      {open && (
        <section aria-label="Optional analytics" className="fixed bottom-4 right-4 left-4 z-50 mx-auto max-w-lg rounded-lg border border-slate-300 bg-white p-4 text-left text-slate-900 shadow-lg">
          <strong>Optional analytics</strong>
          <p className="my-2">Allow Google Analytics to use cookies and send page visits and device information to Google? You can change your choice below the footer.</p>
          <div className="flex gap-3">
            <button type="button" className="rounded border border-slate-400 px-3 py-1" onClick={() => choose("decline")}>Decline</button>
            <button type="button" className="rounded border border-slate-400 px-3 py-1" onClick={() => choose("allow")}>Allow analytics</button>
          </div>
        </section>
      )}
    </div>
  );
}
