## Example File

This directory includes **`HG002Trio.vcf.gz`**, a small subset of variant calls for the **HG002/HG003/HG004 trio** plus **20 unrelated samples** used by `novoCaller`. Variants are fully annotated with **VEP**, and both `novoCaller` and `comHet` were run to refine calls for **HG002**.

This file is intended as a lightweight test dataset for:
- the `granite` **VCF parser** (`granite.lib.vcf_parser`)
- the **filtering utilities** (e.g. `filterByTag`)

---

## Example: VCF parser (child non-reference, parents homozygous reference)

This example extracts variants where the child is **not homozygous reference (0/0)** and **both parents are homozygous reference (0/0)**.

```python
from granite.lib import vcf_parser

# Samples to extract variants for
child = "NA24385_sample"
father = "NA24149_sample"
mother = "NA24143_sample"

# Load the VCF file
vcf_obj = vcf_parser.Vcf("HG002Trio.vcf.gz")

# Open output file and write the header
with open("HG002Trio_test_parser.vcf", "w") as fo:
    vcf_obj.write_header(fo)

    # Iterate through variants
    for vnt_obj in vcf_obj.parse_variants():

        # Extract genotype (GT) values
        GT_child = vnt_obj.get_genotype_value(child, "GT")
        GT_father = vnt_obj.get_genotype_value(father, "GT")
        GT_mother = vnt_obj.get_genotype_value(mother, "GT")

        # Child is non-reference; both parents are homozygous reference
        if GT_child not in ("0/0", "0|0") and GT_father in ("0/0", "0|0") and GT_mother in ("0/0", "0|0"):
            vcf_obj.write_variant(fo, vnt_obj)
```

**Note:**  
The condition above selects variants where the child carries at least one alternate allele (heterozygous or homozygous alternate). If you want to restrict to *only* homozygous alternate calls, explicitly check for `1/1` or `1|1`.

---

## Example: `filterByTag` (high-confidence novoCaller calls)

Extract variants with `novoPP >= 0.9`:

```bash
granite filterByTag \
  -i HG002Trio.vcf.gz \
  -o HG002Trio_test_filter.vcf \
  -t 'novoPP/0.9/>=/float/any'
```

---

## Example: `filterByTag` (novoPP ≥ 0.9 and missense variants)

Extract variants with `novoPP >= 0.9` **and** annotated as `missense_variant` in the VEP `Consequence` field:

```bash
granite filterByTag \
  -i HG002Trio.vcf.gz \
  -o HG002Trio_test_filter.vcf \
  -t 'novoPP/0.9/>=/float/any' \
     'Consequence/missense_variant/~/str/any/field=|/entry=,/value=&' \
  --logic all
```

**Notes:**
- The `~` operator is used to match substrings, which is required for VEP `Consequence` values such as  
  `missense_variant&splice_region_variant`.
- `field=|` specifies the VEP field separator.
- `entry=,` specifies the separator between multiple transcript annotations.
- `value=&` specifies the separator for multiple consequence terms within a transcript.
- All string matching is case-sensitive.
