"use client";

import { useState, useRef, useEffect } from "react";
import { SendHorizonal } from "lucide-react";
import { cn } from "@/lib/utils";
import { requestPlan } from "@/lib/api";
import type { Plan } from "@lattice/sdk/api";

interface Message {
  role: "user" | "assistant";
  content: string;
}

interface ChatPanelProps {
  sessionId: string;
  onPlanUpdate: (plan: Plan) => void;
}

const WELCOME: Message = {
  role: "assistant",
  content:
    'Hello! Upload an .h5ad file, then describe your analysis. For example: "Run standard QC and clustering. This is human PBMC, treated vs. control."',
};

export function ChatPanel({ sessionId, onPlanUpdate }: ChatPanelProps) {
  const [messages, setMessages] = useState<Message[]>([WELCOME]);
  const [input, setInput] = useState("");
  const [loading, setLoading] = useState(false);
  const bottomRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    bottomRef.current?.scrollIntoView({ behavior: "smooth" });
  }, [messages]);

  const send = async () => {
    const text = input.trim();
    if (!text || loading) return;

    const userMsg: Message = { role: "user", content: text };
    setMessages((prev) => [...prev, userMsg]);
    setInput("");
    setLoading(true);

    try {
      const plan = await requestPlan(sessionId, { user_message: text });
      onPlanUpdate(plan);
      const assistantMsg: Message = {
        role: "assistant",
        content: `Generated a ${plan.steps.length}-step plan: ${plan.rationale}`,
      };
      setMessages((prev) => [...prev, assistantMsg]);
    } catch (err) {
      const errMsg: Message = {
        role: "assistant",
        content: `Error: ${err instanceof Error ? err.message : "Request failed."}`,
      };
      setMessages((prev) => [...prev, errMsg]);
    } finally {
      setLoading(false);
    }
  };

  const onKeyDown = (e: React.KeyboardEvent<HTMLTextAreaElement>) => {
    if (e.key === "Enter" && !e.shiftKey) {
      e.preventDefault();
      void send();
    }
  };

  return (
    <div className="flex flex-col h-full" data-testid="chat-panel">
      {/* Messages */}
      <div className="flex-1 overflow-y-auto py-2 space-y-3 pr-1">
        {messages.map((msg, i) => (
          <div
            key={i}
            className={cn(
              "text-xs leading-relaxed rounded px-3 py-2",
              msg.role === "user"
                ? "bg-[hsl(217,91%,60%,0.12)] text-[hsl(210,40%,92%)] ml-4"
                : "text-[hsl(215,20%,70%)]"
            )}
          >
            {msg.role === "assistant" && (
              <span className="text-[10px] uppercase tracking-wider text-[hsl(215,20%,45%)] block mb-1">
                Lattice
              </span>
            )}
            <p className="whitespace-pre-wrap">{msg.content}</p>
          </div>
        ))}
        {loading && (
          <div className="text-xs text-[hsl(215,20%,45%)] px-3 animate-pulse">
            Thinking…
          </div>
        )}
        <div ref={bottomRef} />
      </div>

      {/* Input */}
      <div className="flex items-end gap-2 pt-2 border-t border-[hsl(217,32%,16%)]">
        <textarea
          value={input}
          onChange={(e) => setInput(e.target.value)}
          onKeyDown={onKeyDown}
          placeholder="Describe your analysis…"
          rows={2}
          className={cn(
            "flex-1 resize-none rounded bg-[hsl(217,32%,10%)] border border-[hsl(217,32%,18%)]",
            "px-3 py-2 text-xs text-[hsl(210,40%,92%)] placeholder:text-[hsl(215,20%,40%)]",
            "focus:outline-none focus:border-[hsl(217,91%,60%)] transition-colors"
          )}
          disabled={loading}
          aria-label="Chat input"
        />
        <button
          onClick={() => void send()}
          disabled={loading || !input.trim()}
          className={cn(
            "p-2 rounded bg-[hsl(217,91%,60%)] text-white transition-opacity",
            "disabled:opacity-40 hover:enabled:opacity-90"
          )}
          aria-label="Send message"
        >
          <SendHorizonal className="w-4 h-4" />
        </button>
      </div>
    </div>
  );
}
