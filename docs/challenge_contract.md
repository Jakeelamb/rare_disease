# Challenge operating contract

**Pinned on:** 2026-08-25<br>
**Space revision:** `37e25dceda63ecec7c5b2ebeffd1ea0052ad886e`<br>
**Dataset revision:** `f534cb0c1a607110c6dad0194299bd3dd62df542`

Live drift check on 2026-08-29: Space HEAD was
`1c761cc23d90aebe6a011fd5b0b99517df42408c`. The evaluator, CSV template,
and scoring semantics were unchanged. Track 1 instructions now require users
of a commercial LLM/AI assistant to disclose the provider, plan or tier, and
relevant data-handling setting in the methods report. They also clarify that a
GitHub repository may remain private during the Hackathon but must be public
when it ends.

This is an engineering checklist derived from the pinned official challenge
code and rules. The [full evidence brief](research/challenge_and_methods.md)
contains links and scorer details. Re-check the live rules before every upload.

## Data boundary

- Keep source genomic data and genome-scale intermediate/derived data under
  ignored local paths. Do not grant control rights to unregistered people or
  services.
- Do not recontact the child, family, or MVA Society contacts.
- A hosted service may be used as a processor only when its terms take no
  training/use rights in inputs or outputs and limit retention in time and
  purpose. Record provider, plan/tier, and relevant setting; do not submit
  feedback containing challenge content.
- Delete VCF/BAM/CRAM files, copies/slices/reformats, genotype-scale annotated
  tables, genome-bearing prompts/logs, and genome-trained model artifacts from
  systems under participant control within 30 days after challenge close,
  then send the required deletion confirmation email.
- Ranked candidate variants, HPO terms, gene/pathway rankings, code, reports,
  and leaderboard findings may be retained. Public artifacts must not expose
  a meaningful portion of the child's genome.

## Track 1 release gate

- Exact 12-column CSV schema; `PROBAND01`; GRCh38; at most 10 rows.
- Emit `chrN` chromosome labels even though the source VCF uses `N`.
- Put an exact proposed compound-heterozygous pair in one row.
- Keep EPCR values unique and strictly descending; review every row manually.
- Do not place incidental findings above primary hypotheses because the pinned
  evaluator does not exclude `finding_type` from automated scoring.
- Attach the GitHub repository and a PDF/Markdown methods report only after a
  privacy review. The repository may remain private during the Hackathon but
  must be public when it ends.
- If an LLM/AI assistant was used, record its provider, plan/tier, and relevant
  data-handling setting in the methods report.

## Track 2 release gate

- Establish the genetic and molecular mechanism before nominating a medicine.
- Restrict candidates to market-approved medicines and distinguish mechanistic
  plausibility from efficacy, safety, or clinical advice.
- Provide claim-level source identifiers, release dates, and contradictory
  evidence; include a validation plan and uncertainty statement.
- Submit one report, one public repository, and one three-minute pitch through
  the designated team member.

## Publication and license

Participant code, reports, predictions, and results are CC BY 4.0 under the
challenge rules. The controlled patient data are not redistributable. Observe
the post-close publication embargo until organizers publicly end it.
