#!/usr/bin/env bash
# Reproduction driver.
#
#   bash run_experiments.sh          verify the three exact-rational
#                                    certificates against their frozen
#                                    outputs and regenerate the
#                                    derived-numbers file (minutes)
#   bash run_experiments.sh --full   additionally rerun every
#                                    experiment from the raw data
#                                    (hours)
#
# Paths resolve relative to this repository; override with COLAS_BASE.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
export COLAS_BASE="$HERE"
cd "$HERE"

echo "== [1/2] certificates (exact rational; frozen-output comparison) =="
python3 code/colas_certificate_v2.py | tee /tmp/cert_v2.out
diff <(grep -v '^#' /tmp/cert_v2.out) \
     <(grep -v '^#' code/cert_v2_expected_output.txt) \
  && echo "certificate v2: output matches frozen"
python3 code/colas_certificate_d2.py | tee /tmp/cert_d2.out
diff <(grep -v '^#' /tmp/cert_d2.out) \
     <(grep -v '^#' code/cert_d2_expected_output.txt) \
  && echo "certificate d2: output matches frozen"
python3 code/ieff_face_certificate.py | tee /tmp/cert_ieff.out
diff /tmp/cert_ieff.out code/ieff_cert_expected_output.txt \
  && echo "certificate I_eff face: output matches frozen"

if [[ "${1:-}" == "--full" ]]; then
  echo "== [full] simulations =="
  python3 code/run_experiments.py
  python3 code/sim_extended.py
  python3 code/sim_spec_test.py
  echo "== [full] eleven-city pipeline (both platforms) =="
  python3 code/brightkite_cities.py
  python3 code/gowalla_cities.py
  python3 code/cities_gof.py
  python3 code/cities_distnull.py
  python3 code/sf_extras.py
  python3 code/design_functionals.py
  python3 code/ieff_profile.py
  echo "== [full] augmented class and competing models =="
  python3 code/rich_fit_city.py "SF Bay" brightkite
  RICH_FIXB=dir python3 code/rich_fit_city.py "Austin" gowalla
  python3 code/competitors_sf.py
  echo "== [full] collaboration network (OpenAlex) =="
  OA_BASE=1 python3 code/oa_analysis.py
  OA_SENS=1 python3 code/oa_analysis.py
  python3 code/oa_rich.py free
  python3 code/oa_rich.py fixdir
  python3 code/oa_rich.py free2 '[1.0,0.8,0.3,40.0,1.5,0.02,0.1]'
fi

echo "== [2/2] derived numbers from frozen results =="
python3 code/make_numbers.py
echo "done: revision/numbers.tex regenerated from runs/"
