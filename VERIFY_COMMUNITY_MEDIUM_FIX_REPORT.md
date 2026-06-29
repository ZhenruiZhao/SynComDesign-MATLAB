# Verify Community Medium Fix Report

Project checked: `D:\SynComDesign-MATLAB-community-medium-fix`

Final classification: **A. external-medium-only restriction, cross-feeding allowed**

Final status: **COMMUNITY_MEDIUM_FIX_VERIFIED**

## Findings

- External shared exchange present: YES
- Major medium entries map to shared exchange: YES
- Unlisted external uptake closed: YES
- Strain-interface changed by medium: NO
- Internal transport changed by medium: NO
- Cross-feeding structurally allowed: YES
- Default mode external_shared_only: YES
- Transport reactions misclassified as interface: NO

## Medium Application Order

The default path builds the community, creates shared external exchanges, then calls `applyCommunityExternalMedium`. Old `applyMedium` is used only for explicit `legacy_all_exchange` community mode or temporary single-strain validation.

## Reaction Classification

The verification table includes all audited reactions. Transport reactions such as `NO3tex`, `NO2tex`, `O2tex`, and `N2Otpp` are classified as `internal_transport` with `medium_applies = false` in the audited community models.

## Medium Mapping and Bounds

Medium mapping found fraction: 0.96
Unlisted external exchanges were checked for closed uptake. Interface and transport before/after tables verify whether medium changed their bounds.

## Cross-feeding

The structural audit checks non-medium compounds for closed external uptake plus open producer/consumer strain interfaces and a shared metabolite.

## Biomass

Legacy vs fixed growth comparison rows: 31
Zero biomass in fixed mode is acceptable only when not caused by medium changes to interface or transport reactions.

## Python Alignment Target

Can be used as Python alignment target: YES

## Remaining Non-Blocking Warnings

The only unresolved medium mappings in the audited combinations are `R_EX_ribflv_e` and `R_EX_slnt_e`; their expected shared exchanges are absent from those community models. This is reported explicitly in `10_medium_warnings.tsv` and is not the old broad `Medium exchange reaction not found: R_EX_*` failure pattern.

## Issues Still Requiring Attention

No blocking issue found by this verification.

## Manual Checks

Review remaining missing medium mappings and decide whether the underlying strains are expected to contain those shared metabolites.
