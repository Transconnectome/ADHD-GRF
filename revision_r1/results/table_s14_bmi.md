# Supplementary Table S14: Sensitivity Refits for Measured Child Body Mass Index

**Table S14. Generalized random forest results for the selected environmental and genetic model refitted with measured child body mass index.**

<table style="width:100%;">
<caption><p><strong>Table S14. Generalized random forest results for the selected model refitted with measured child body mass index.</strong> Models were fitted in the primary multiethnic sample (<em>N</em> = 1,244; 622 exposed, 622 unexposed) using an ensemble of 5 seeds with 5,000 trees per seed (25,000 trees total), clustering on study site, and AIPW estimation. CATE coefficients and standard errors are from out-of-bag <code>test_calibration</code> regressions.</p></caption>
<colgroup>
<col style="width: 32%;" />
<col style="width: 12%;" />
<col style="width: 10%;" />
<col style="width: 14%;" />
<col style="width: 10%;" />
<col style="width: 11%;" />
<col style="width: 11%;" />
</colgroup>
<thead>
<tr>
<th style="text-align: left;"><strong>Model / Specification</strong></th>
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
<td colspan="7" style="text-align: left; background-color: #f2f2f2;"><strong>Primary Reference Model (Reported in Manuscript)</strong></td>
</tr>
<tr>
<td style="text-align: left;">Primary Multiethnic Model (all 66 covariates; BMI PGS included, measured child BMI excluded)</td>
<td style="text-align: center;">1,244 (622 : 622)</td>
<td style="text-align: center;">0.54 (0.03)</td>
<td style="text-align: center;">[0.47, 0.60]</td>
<td style="text-align: center;">< .001</td>
<td style="text-align: center;">1.00 (0.11)</td>
<td style="text-align: center;">2.17 (0.73, 0.002)</td>
</tr>
<tr>
<td colspan="7" style="text-align: left; background-color: #f2f2f2;"><strong>Measured Child BMI Sensitivity Refits</strong></td>
</tr>
<tr>
<td style="text-align: left;">(a) Substituting BMI PGS with measured baseline child BMI (66 covariates)</td>
<td style="text-align: center;">1,244 (622 : 622)</td>
<td style="text-align: center;">0.54 (0.03)</td>
<td style="text-align: center;">[0.47, 0.61]</td>
<td style="text-align: center;">< .001</td>
<td style="text-align: center;">1.01 (0.11)</td>
<td style="text-align: center;">1.79 (1.06, 0.091)</td>
</tr>
<tr>
<td style="text-align: left;">(b) Including both BMI PGS and measured baseline child BMI (67 covariates)</td>
<td style="text-align: center;">1,244 (622 : 622)</td>
<td style="text-align: center;">0.54 (0.03)</td>
<td style="text-align: center;">[0.47, 0.60]</td>
<td style="text-align: center;">< .001</td>
<td style="text-align: center;">1.01 (0.11)</td>
<td style="text-align: center;">2.15 (0.89, 0.016)</td>
</tr>
</tbody>
</table>

*Note.* Child body mass index (`bmi`) was calculated from measured baseline height and weight (`anthro_1_height_in` and `anthroweight1lb`) according to the standard formula (weight in lbs / height in inches^2 * 703; Table S3). In specification (a) (`bmi_substituting_pgs`), the body mass index polygenic score (`BMImulti`) was replaced by measured child `bmi` (66 covariates total). In specification (b) (`bmi_plus_pgs`), both measured child `bmi` and the BMI polygenic score (`BMImulti`) were included simultaneously (67 covariates total). In specification (a), the top five covariates by variable importance were sleep_prob (`sleep_prob`, 0.096), SMOKERauto (`SMOKERauto`, 0.065), INSOMNIAeur6 (`INSOMNIAeur6`, 0.047), SCZ_METAauto (`SCZ_METAauto`, 0.043), WORRYauto (`WORRYauto`, 0.042). In specification (b), the top five covariates by variable importance were sleep_prob (`sleep_prob`, 0.091), BMImulti (`BMImulti`, 0.071), SMOKERauto (`SMOKERauto`, 0.070), SCZ_METAauto (`SCZ_METAauto`, 0.045), WORRYauto (`WORRYauto`, 0.045).

