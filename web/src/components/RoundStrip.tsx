import { useIntegration } from '../integrationStore'
import { roundStripText } from '../mergeRounds'

/** Section 2: the segmentation base (original → rounds) with undo and back-to-original. */
export default function RoundStrip() {
  const ig = useIntegration()
  return <div className="round-strip">
    <p className="small" aria-label="Segmentation base"><b>Segmentation base:</b>{' '}
      {['Original', ...ig.rounds.map((r, i) => `Round ${i + 1} (${r.label}, ${r.merges} merges)`)].map((t, i, a) =>
        <span key={i} title={i ? JSON.stringify(ig.rounds[i - 1].rules, null, 1) : 'original accepted breakpoints'} className={i === a.length - 1 ? 'round-current' : ''}>{i ? ' → ' : ''}{t}</span>)}
    </p>
    <span className="sr-only">{roundStripText(ig.rounds)}</span>
    <button className="btn-sm" disabled={!ig.rounds.length} onClick={() => ig.undoRound()}>Undo last round</button>
    <button className="btn-sm" disabled={!ig.rounds.length} onClick={() => ig.resetRounds()}>Back to original</button>
    {ig.roundsNotice && <p className="notice small" role="status">{ig.roundsNotice}</p>}
    {ig.baseError && <p className="notice small" role="status">{ig.baseError}</p>}
  </div>
}
