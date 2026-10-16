# Supplementary Table S16: Causal Forest Sensitivity Refitting in the European Sample Without Oversampling

**Table S16. Selected model (environmental and genetic covariates) refitted in the European-ancestry sample before oversampling, under two definitions of the unexposed group.**

<table style="width:100%;">
<caption><p><strong>Table S16. Causal forest results in the non-oversampled European sample compared to the oversampled European analysis and primary multiethnic analysis.</strong> Models were fitted using an ensemble of 5 seeds with 5,000 trees per seed (25,000 trees total), clustering on study site, and AIPW estimation. CATE coefficients and standard errors are from out-of-bag <code>test_calibration</code> regressions.</p></caption>
<colgroup>
<col style="width: 28%;" />
<col style="width: 10%;" />
<col style="width: 10%;" />
<col style="width: 14%;" />
<col style="width: 10%;" />
<col style="width: 14%;" />
<col style="width: 14%;" />
</colgroup>
<thead>
<tr>
<th style="text-align: left;"><strong>Model / Sample Definition</strong></th>
<th style="text-align: center;"><strong><em>N</em> (Treated : Control)</strong></th>
<th style="text-align: center;"><strong>ATE (SE)</strong></th>
<th style="text-align: center;"><strong>ATE 95% CI</strong></th>
<th style="text-align: center;"><strong>ATE <em>p</em></strong></th>
<th style="text-align: center;"><strong>CATE Mean Pred. (SE)</strong></th>
<th style="text-align: center;"><strong>CATE Diff. Slope (SE, <em>p</em>)</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td colspan="7" style="text-align: left; background-color: #f2f2f2;"><strong>Primary Reference Analyses (Reported in Manuscript)</strong></td>
</tr>
<tr>
<td style="text-align: left;">Primary Multiethnic Sample (no synthetic observations)</td>
<td style="text-align: center;">1,244 (622 : 622)</td>
<td style="text-align: center;">0.54 (0.03)</td>
<td style="text-align: center;">[0.47, 0.60]</td>
<td style="text-align: center;">< .001</td>
<td style="text-align: center;">1.00 (0.11)</td>
<td style="text-align: center;">2.17 (0.73, 0.002)</td>
</tr>
<tr>
<td style="text-align: left;">European post-SMOTE (secondary; contains 257 synthetic treated)</td>
<td style="text-align: center;">1,239 (620 : 619)</td>
<td style="text-align: center;">0.52 (0.05)</td>
<td style="text-align: center;">[0.41, 0.63]</td>
<td style="text-align: center;">< .001</td>
<td style="text-align: center;">1.05 (0.17)</td>
<td style="text-align: center;">2.55 (0.67, < .001)</td>
</tr>
<tr>
<td colspan="7" style="text-align: left; background-color: #f2f2f2;"><strong>P12 (R3-5): Non-Oversampled European Sensitivity Refits (Original 363 Exposed Children)</strong></td>
</tr>
<tr>
<td style="text-align: left;">(a) Sex-matched unexposed group (downsampled 206 M : 157 F)</td>
<td style="text-align: center;">726 (363 : 363)</td>
<td style="text-align: center;">0.57 (0.05)</td>
<td style="text-align: center;">[0.48, 0.66]</td>
<td style="text-align: center;">< .001</td>
<td style="text-align: center;">0.97 (0.12)</td>
<td style="text-align: center;">1.84 (0.93, 0.047)</td>
</tr>
<tr>
<td style="text-align: left;">(b) Full unexposed group (all available complete controls)</td>
<td style="text-align: center;">6,106 (363 : 5,743)</td>
<td style="text-align: center;">0.53 (0.06)</td>
<td style="text-align: center;">[0.42, 0.64]</td>
<td style="text-align: center;">< .001</td>
<td style="text-align: center;">1.10 (0.14)</td>
<td style="text-align: center;">1.10 (0.54, 0.041)</td>
</tr>
</tbody>
</table>

### Key Methodological Notes:
1. **P12(a) Sex-Matched Specification**: Comprises the original 363 exposed children of European ancestry together with a sex-matched draw of 363 unexposed children (matching the 206 male / 157 female ratio of the exposed group, random seed 1), total *N* = 726. In this sample, the differential forest prediction slope was 1.8439 (SE = 0.9253, two-sided *p* = 0.047), which is nominally significant at alpha = 0.05 but below the minimum detectable effect threshold (MDE = 2.5908 at 80% power).
2. **P12(b) Full Control Specification**: Comprises the original 363 exposed children together with all 5,743 available unexposed European children (*N* = 6,106). In this full sample, the average treatment effect is 0.5298 (SE = 0.0557, *p* < .001), and differential prediction slope is 1.1038 (SE = 0.5395, two-sided *p* = 0.041).
3. **Data Governance & Reproducibility**: All models were fitted using R package `grf` (version 2.6.1) with 25,000 trees across 5 independent seeds with study site clustering. No individual-level ABCD data are exposed.

