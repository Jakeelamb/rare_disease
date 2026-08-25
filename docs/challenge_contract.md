# Challenge operating contract

**Pinned on:** 2026-08-25<br>
**Space revision:** `37e25dceda63ecec7c5b2ebeffd1ea0052ad886e`<br>
**Dataset revision:** `f534cb0c1a607110c6dad0194299bd3dd62df542`

This is an engineering checklist derived from the pinned official challenge
code and rules. The [full evidence brief](research/challenge_and_methods.md)
contains links and scorer details. Re-check the live rules before every upload.

## Data boundary

- Keep all source, intermediate, and derived patient data under ignored local
  paths. Do not grant access to unregistered people or hosted services.
- Do not recontact the child, family, or MVA Society contacts.
- Treat patient-level calls, HPO rows, candidate lists, and model prompts as
  restricted. Hosted API use is paused pending organizer clarification.
- Delete source and derived data from every environment within 30 days after
  challenge close, then send the required deletion confirmation email.
- Public artifacts may contain code, configuration, aggregate QC, cited public
  knowledge, and carefully minimized conclusions only.

## Track 1 release gate

- Exact 12-column CSV schema; `PROBAND01`; GRCh38; at most 10 rows.
- Emit `chrN` chromosome labels even though the source VCF uses `N`.
- Put an exact proposed compound-heterozygous pair in one row.
- Keep EPCR values unique and strictly descending; review every row manually.
- Do not place incidental findings above primary hypotheses because the pinned
  evaluator does not exclude `finding_type` from automated scoring.
- Attach the public code repository and a methods report only after a privacy
  review confirms neither contains restricted data.

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
