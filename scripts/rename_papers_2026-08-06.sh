#!/usr/bin/env bash
# GENERATED 2026-08-06 by the library-curation pass. One-off.
#
# Renames every PDF in ~/papers/pymrm-gallery/ to the repo scheme
#   papers: <FirstAuthor><Year>-<slug>-<JournalAbbrev><Vol>-<FirstPage>.pdf
#   books:  BOOK-<Authors>-<Year>-<short-title>-<edition>.pdf
#   non-case files kept for the record: MISC-<...>.pdf
#
# Every identity behind a name is recorded in docs/papers-inventory.yaml,
# together with the old name, so provenance chains survive.
#
# Run scripts/update_paper_references_2026-08-06.py AFTER this.
set -euo pipefail
P="${PYMRM_PAPERS:-$HOME/papers/pymrm-gallery}"
cd "$P"

# --- byte-identical duplicates: moved aside, not deleted ---------------------
mkdir -p duplicates
[ -f "Biotech   Bioengineering - November 1968 - Andrews.pdf" ] && mv -n -- "Biotech   Bioengineering - November 1968 - Andrews.pdf" duplicates/ || true
[ -f "Journal of Biochemical and Microbiological Technology and Engineering - December 1959 - Luedeking - A kinetic study of the.pdf" ] && mv -n -- "Journal of Biochemical and Microbiological Technology and Engineering - December 1959 - Luedeking - A kinetic study of the.pdf" duplicates/ || true
[ -f "Journal of Biochemical and Microbiological Technology and Engineering - Luedeking.pdf" ] && mv -n -- "Journal of Biochemical and Microbiological Technology and Engineering - Luedeking.pdf" duplicates/ || true

# --- renames ----------------------------------------------------------------
[ -f "ef070025k.pdf" ] && mv -n -- "ef070025k.pdf" "Abad2007-oxygen-carrier-reduction-kinetics-EnergyFuels21-1843.pdf" || true
[ -f "Biotech   Bioengineering - November 1968 - Andrews - A mathematical model for the continuous culture of microorganisms.pdf" ] && mv -n -- "Biotech   Bioengineering - November 1968 - Andrews - A mathematical model for the continuous culture of microorganisms.pdf" "Andrews1968-substrate-inhibition-continuous-culture-BiotechnolBioeng10-707.pdf" || true
[ -f "1-s2.0-0009250957850283-main.pdf" ] && mv -n -- "1-s2.0-0009250957850283-main.pdf" "Aris1957-shape-generalised-modulus-CES6-262-PREVIEW1P.pdf" || true
[ -f "Bird R B, Stewart W E, Lightfoot E N (2nd Ed,Wiley,2002) - Transport Phenomena.pdf" ] && mv -n -- "Bird R B, Stewart W E, Lightfoot E N (2nd Ed,Wiley,2002) - Transport Phenomena.pdf" "BOOK-BirdStewartLightfoot-2002-transport-phenomena-2ed.pdf" || true
[ -f "cambridge-mathematical-library-sydney-chapman-t.-g.-cowling-c.-cercignani-the-mathematical-theory-of-non-uniform-gases_-an-account-of-the-kinetic-theory-of-viscosity-thermal-conduction-a.pdf" ] && mv -n -- "cambridge-mathematical-library-sydney-chapman-t.-g.-cowling-c.-cercignani-the-mathematical-theory-of-non-uniform-gases_-an-account-of-the-kinetic-theory-of-viscosity-thermal-conduction-a.pdf" "BOOK-ChapmanCowling-1970-non-uniform-gases-3ed.pdf" || true
[ -f "Froment_Bischoff_Chemical_Reactor_Analysis_and_Design.pdf" ] && mv -n -- "Froment_Bischoff_Chemical_Reactor_Analysis_and_Design.pdf" "BOOK-FromentDeWildeBischoff-2011-chemical-reactor-analysis-and-design-3ed.pdf" || true
[ -f "2015.205681.Industrial-Chemical.pdf" ] && mv -n -- "2015.205681.Industrial-Chemical.pdf" "BOOK-HougenWatson-1936-industrial-chemical-calculations-2ed.pdf" || true
[ -f "process-calculation-by-watson.pdf" ] && mv -n -- "process-calculation-by-watson.pdf" "BOOK-HougenWatson-1947-chemical-process-principles-combined.pdf" || true
[ -f "Chemical_reaction_engineering-levenspiel.pdf" ] && mv -n -- "Chemical_reaction_engineering-levenspiel.pdf" "BOOK-Levenspiel-1999-chemical-reaction-engineering-3ed.pdf" || true
[ -f "7981631866938162.pdf" ] && mv -n -- "7981631866938162.pdf" "BOOK-LiKwauk-1994-particle-fluid-two-phase-flow-EMMS-1ed.pdf" || true
[ -f "CHEMREACFUN-book-2nd-edition-5th-printing.pdf" ] && mv -n -- "CHEMREACFUN-book-2nd-edition-5th-printing.pdf" "BOOK-RawlingsEkerdt-2025-chemical-reactor-analysis-design-fundamentals-2ed-5pr.pdf" || true
[ -f "RossTaylorR.Krishna-MulticomponentMassTransferWileySeriesinChemicalEngineering1993.pdf" ] && mv -n -- "RossTaylorR.Krishna-MulticomponentMassTransferWileySeriesinChemicalEngineering1993.pdf" "BOOK-TaylorKrishna-1993-multicomponent-mass-transfer-1ed.pdf" || true
[ -f "1-s2.0-0300946789850026-main.pdf" ] && mv -n -- "1-s2.0-0300946789850026-main.pdf" "Baldyga1989-micromixing-simplification-pt1-ChemEngJ42-83.pdf" || true
[ -f "1-s2.0-0300946789850038-main.pdf" ] && mv -n -- "1-s2.0-0300946789850038-main.pdf" "Baldyga1989-micromixing-simplification-pt2-ChemEngJ42-93.pdf" || true
[ -f "i160071a009.pdf" ] && mv -n -- "i160071a009.pdf" "Beeckman1979-site-coverage-pore-blockage-IECFund18-245.pdf" || true
[ -f "1-s2.0-S0098135403002072-main.pdf" ] && mv -n -- "1-s2.0-S0098135403002072-main.pdf" "Bezzo2004-hybrid-multizonal-cfd-pt2-CACE28-513.pdf" || true
[ -f "AIChE Journal - May 1980.pdf" ] && mv -n -- "AIChE Journal - May 1980.pdf" "Bhatia1980-random-pore-model-pt1-AIChEJ26-379.pdf" || true
[ -f "Literature-PredictionofMassTransferColumns.pdf" ] && mv -n -- "Literature-PredictionofMassTransferColumns.pdf" "Billet1999-mass-transfer-columns-packings-TransIChemE77A-498.pdf" || true
[ -f "AIChE Journal - March 1965 - Bischoff - Effectiveness factors for general reaction rate forms.pdf" ] && mv -n -- "AIChE Journal - March 1965 - Bischoff - Effectiveness factors for general reaction rate forms.pdf" "Bischoff1965-generalised-effectiveness-factor-AIChEJ11-351.pdf" || true
[ -f "F3_5.pdf" ] && mv -n -- "F3_5.pdf" "Bosch1989-CO2-absorption-promoted-carbonate-CES44-2735.pdf" || true
[ -f "ja01269a023.pdf" ] && mv -n -- "ja01269a023.pdf" "Brunauer1938-BET-multimolecular-adsorption-JACS60-309.pdf" || true
[ -f "1-s2.0-000925096187005X-main.pdf" ] && mv -n -- "1-s2.0-000925096187005X-main.pdf" "Calderbank1961-mass-transfer-agitated-vessels-CES16-39.pdf" || true
[ -f "AIChE Journal - March 1960 - Charlesworth - Evaporation from drops containing dissolved solids.pdf" ] && mv -n -- "AIChE Journal - March 1960 - Charlesworth - Evaporation from drops containing dissolved solids.pdf" "Charlesworth1960-evaporation-drops-dissolved-solids-AIChEJ6-9.pdf" || true
[ -f "ma00237a002.pdf" ] && mv -n -- "ma00237a002.pdf" "Chiu1983-gel-effect-free-radical-polymerization-Macromolecules16-348.pdf" || true
[ -f "ie8b02111.pdf" ] && mv -n -- "ie8b02111.pdf" "Criado2018-CaO-carbonation-temperature-IECR57-12595.pdf" || true
[ -f "ie50498a055.pdf" ] && mv -n -- "ie50498a055.pdf" "Danckwerts1951-liquid-film-coefficients-gas-absorption-IEC43-1460.pdf" || true
[ -f "1-s2.0-0009250953800011-main.pdf" ] && mv -n -- "1-s2.0-0009250953800011-main.pdf" "Danckwerts1953-continuous-flow-systems-CES2-1.pdf" || true
[ -f "1-s2.0-S0360128507000214-main.pdf" ] && mv -n -- "1-s2.0-S0360128507000214-main.pdf" "DiBlasi2008-wood-biomass-pyrolysis-review-PECS34-47.pdf" || true
[ -f "AIChE Journal - July 1979 - Dixon - Theoretical prediction of effective heat transfer parameters in packed beds.pdf" ] && mv -n -- "AIChE Journal - July 1979 - Dixon - Theoretical prediction of effective heat transfer parameters in packed beds.pdf" "Dixon1979-effective-heat-transfer-parameters-AIChEJ25-663.pdf" || true
[ -f "Doyle_1993_J._Electrochem._Soc._140_1526.pdf" ] && mv -n -- "Doyle_1993_J._Electrochem._Soc._140_1526.pdf" "Doyle1993-lithium-cell-galvanostatic-model-JES140-1526.pdf" || true
[ -f "AIChE Journal - March 1962 - Duncan - An experimental study of three component gas diffusion.pdf" ] && mv -n -- "AIChE Journal - March 1962 - Duncan - An experimental study of three component gas diffusion.pdf" "Duncan1962-three-component-gas-diffusion-AIChEJ8-38.pdf" || true
[ -f "i160028a013.pdf" ] && mv -n -- "i160028a013.pdf" "Dyson1968-ammonia-synthesis-kinetics-diffusion-IECFund7-605.pdf" || true
[ -f "1-s2.0-0009250968870563-main.pdf" ] && mv -n -- "1-s2.0-0009250968870563-main.pdf" "Edwards1968-gas-dispersion-packed-beds-CES23-109.pdf" || true
[ -f "211103107pdfcreator.pdf" ] && mv -n -- "211103107pdfcreator.pdf" "Ergun1952-fluid-flow-through-packed-columns-ChemEngProg48-89.pdf" || true
[ -f "i160046a001.pdf" ] && mv -n -- "i160046a001.pdf" "Feng1973-isothermal-diffusion-porous-solids-IECFund12-143.pdf" || true
[ -f "ef00034a011.pdf" ] && mv -n -- "ef00034a011.pdf" "Fletcher1992-CPD-model-pt3-13C-NMR-EnergyFuels6-414.pdf" || true
[ -f "1-s2.0-0009250964850922-main.pdf" ] && mv -n -- "1-s2.0-0009250964850922-main.pdf" "Franckaerts1964-ethanol-dehydrogenation-kinetics-CES19-807.pdf" || true
[ -f "1-s2.0-0009250961800304-main.pdf" ] && mv -n -- "1-s2.0-0009250961800304-main.pdf" "Froment1961-fixed-bed-fouling-pt1-CES16-189.pdf" || true
[ -f "ie50686a006.pdf" ] && mv -n -- "ie50686a006.pdf" "Froment1967-fixed-bed-reactors-design-status-IEC59-18.pdf" || true
[ -f "ie50677a007.pdf" ] && mv -n -- "ie50677a007.pdf" "Fuller1966-diffusion-volumes-IEC58-18.pdf" || true
[ -f "1-s2.0-0032591073800373-main.pdf" ] && mv -n -- "1-s2.0-0032591073800373-main.pdf" "Geldart1973-types-of-gas-fluidization-PowderTechnol7-285.pdf" || true
[ -f "AIChE Journal - 2004 - Gheorghiu - Optimal bimodal pore networks for heterogeneous catalysis.pdf" ] && mv -n -- "AIChE Journal - 2004 - Gheorghiu - Optimal bimodal pore networks for heterogeneous catalysis.pdf" "Gheorghiu2004-optimal-bimodal-pore-networks-AIChEJ50-812.pdf" || true
[ -f "tf9555101540.pdf" ] && mv -n -- "tf9555101540.pdf" "Glueckauf1955-linear-driving-force-TransFaradaySoc51-1540.pdf" || true
[ -f "1-s2.0-0009250988851273-main.pdf" ] && mv -n -- "1-s2.0-0009250988851273-main.pdf" "Graaf1988-methanol-synthesis-kinetics-CES43-3185.pdf" || true
[ -f "Annalen der Physik - 1882 - Graetz - Ueber die W rmeleitungsf higkeit von Fl ssigkeiten.pdf" ] && mv -n -- "Annalen der Physik - 1882 - Graetz - Ueber die W rmeleitungsf higkeit von Fl ssigkeiten.pdf" "Graetz1882-waermeleitungsfaehigkeit-fluessigkeiten-AnnPhys254-79.pdf" || true
[ -f "ef00014a011.pdf" ] && mv -n -- "ef00014a011.pdf" "Grant1989-chemical-percolation-devolatilization-CPD-EnergyFuels3-175.pdf" || true
[ -f "1-s2.0-0009250987850662-main.pdf" ] && mv -n -- "1-s2.0-0009250987850662-main.pdf" "Gunn1987-axial-radial-dispersion-fixed-beds-CES42-363.pdf" || true
[ -f "1-s2.0-025527019380020H-main.pdf" ] && mv -n -- "1-s2.0-025527019380020H-main.pdf" "Gunn1993-axial-dispersion-note-ChemEngProcess32-333.pdf" || true
[ -f "AIChE Journal - May 1976 - Heck.pdf" ] && mv -n -- "AIChE Journal - May 1976 - Heck.pdf" "Heck1976-monolithic-catalyst-modeling-AIChEJ22-477.pdf" || true
[ -f "1-s2.0-0043135487900583-main.pdf" ] && mv -n -- "1-s2.0-0043135487900583-main.pdf" "Henze1987-activated-sludge-model-1-WaterRes21-505.pdf" || true
[ -f "AIChE Journal - 1983 - Herskowitz.pdf" ] && mv -n -- "AIChE Journal - 1983 - Herskowitz.pdf" "Herskowitz1983-trickle-bed-partial-wetting-AIChEJ29-1.pdf" || true
[ -f "1-s2.0-0009250964850478-main.pdf" ] && mv -n -- "1-s2.0-0009250964850478-main.pdf" "Hulburt1964-population-balance-CES19-555.pdf" || true
[ -f "AIChE Journal - September 1987 - Itoh - A membrane reactor using palladium.pdf" ] && mv -n -- "AIChE Journal - September 1987 - Itoh - A membrane reactor using palladium.pdf" "Itoh1987-palladium-membrane-reactor-AIChEJ33-1576.pdf" || true
[ -f "20100036467.pdf" ] && mv -n -- "20100036467.pdf" "Kandula2010-effective-thermal-conductivity-packed-beds-NASA-KSC-2010-007.pdf" || true
[ -f "90131.pdf" ] && mv -n -- "90131.pdf" "Kiani2024-pair-sites-langmuir-hinshelwood-ACSCatal14-10260.pdf" || true
[ -f "ac60131a045.pdf" ] && mv -n -- "ac60131a045.pdf" "Kissinger1957-reaction-kinetics-in-DTA-AnalChem29-1702.pdf" || true
[ -f "1-s2.0-S008207847780341X-main.pdf" ] && mv -n -- "1-s2.0-S008207847780341X-main.pdf" "Kobayashi1977-coal-devolatilization-high-T-SympCombust16-411.pdf" || true
[ -f "AIChE_Journal-1996-Krishna.pdf" ] && mv -n -- "AIChE_Journal-1996-Krishna.pdf" "Krishna1996-bubble-column-gas-holdup-AIChEJ42-2627.pdf" || true
[ -f "1-s2.0-S0009250996004587-main.pdf" ] && mv -n -- "1-s2.0-S0009250996004587-main.pdf" "Krishna1997-maxwell-stefan-review-CES52-861.pdf" || true
[ -f "i260028a001.pdf" ] && mv -n -- "i260028a001.pdf" "Kunii1968-bubbling-bed-model-IECFund7-481.pdf" || true
[ -f "ja02242a004.pdf" ] && mv -n -- "ja02242a004.pdf" "Langmuir1918-adsorption-plane-surfaces-JACS40-1361.pdf" || true
[ -f "AIChE Journal - June 1961 - Larkins.pdf" ] && mv -n -- "AIChE Journal - June 1961 - Larkins.pdf" "Larkins1961-two-phase-concurrent-flow-packed-beds-AIChEJ7-231.pdf" || true
[ -f "AIChE Journal - July 1984 - Lehrer.pdf" ] && mv -n -- "AIChE Journal - July 1984 - Lehrer.pdf" "Lehrer1984-turbulent-axial-dispersion-bubble-column-AIChEJ30-654.pdf" || true
[ -f "1-s2.0-0021951772902278-main.pdf" ] && mv -n -- "1-s2.0-0021951772902278-main.pdf" "Levenspiel1972-deactivating-catalyst-rate-equation-JCatal25-265.pdf" || true
[ -f "Journal of Biochem and Microbiol Techn and Engin - Luedeking.pdf" ] && mv -n -- "Journal of Biochem and Microbiol Techn and Engin - Luedeking.pdf" "Luedeking1959-lactic-acid-fermentation-kinetics-JBMTE1-393.pdf" || true
[ -f "2012083364.pdf" ] && mv -n -- "2012083364.pdf" "MISC-ChemEngProgress-2012-08-shale-gas-issue-CEP108-8.pdf" || true
[ -f "A1993LB49000001.pdf" ] && mv -n -- "A1993LB49000001.pdf" "MISC-Ranz1993-citation-classic-commentary-CurrContents22.pdf" || true
[ -f "ie50546a056.pdf" ] && mv -n -- "ie50546a056.pdf" "MISC-Wilke1955-binary-diffusion-coefficient-estimation-IEC47-1253.pdf" || true
[ -f "1-s2.0-S0920586199000826-main.pdf" ] && mv -n -- "1-s2.0-S0920586199000826-main.pdf" "Maretto1999-slurry-bubble-column-FT-CatalToday52-279.pdf" || true
[ -f "413a375.pdf" ] && mv -n -- "413a375.pdf" "Markos1987-catalyst-deactivation-parameter-estimation-pt4-ChemPap41-375.pdf" || true
[ -f "Marquis_2019_J._Electrochem._Soc._166_A3693.pdf" ] && mv -n -- "Marquis_2019_J._Electrochem._Soc._166_A3693.pdf" "Marquis2019-single-particle-model-electrolyte-JES166-A3693.pdf" || true
[ -f "1-s2.0-S0009250954800054-main.pdf" ] && mv -n -- "1-s2.0-S0009250954800054-main.pdf" "Mars1954-vanadium-oxide-oxidation-CESSuppl3-41.pdf" || true
[ -f "Chemie Ingenieur Technik - December 1993 - Martin - Radiale W rmeleitung in durchstr mten Sch ttungsrohren.pdf" ] && mv -n -- "Chemie Ingenieur Technik - December 1993 - Martin - Radiale W rmeleitung in durchstr mten Sch ttungsrohren.pdf" "Martin1993-radiale-waermeleitung-schuettungsrohren-CIT65-1468.pdf" || true
[ -f "i260040a020.pdf" ] && mv -n -- "i260040a020.pdf" "Mears1971-tests-for-transport-limitations-IECPDD10-541.pdf" || true
[ -f "1-s2.0-S0009250912007099-main.pdf" ] && mv -n -- "1-s2.0-S0009250912007099-main.pdf" "Mills2013-two-dimensional-stefan-tube-CES90-130.pdf" || true
[ -f "1-s2.0-S0032591024012488-main.pdf" ] && mv -n -- "1-s2.0-S0032591024012488-main.pdf" "Mou2025-simplified-ZBS-conductivity-PowderTechnol453-120604.pdf" || true
[ -f "AIChE Journal - January 1965 - Myers.pdf" ] && mv -n -- "AIChE Journal - January 1965 - Myers.pdf" "Myers1965-thermodynamics-mixed-gas-adsorption-IAST-AIChEJ11-121.pdf" || true
[ -f "1-s2.0-0022024868901796-main.pdf" ] && mv -n -- "1-s2.0-0022024868901796-main.pdf" "Nyvlt1968-nucleation-kinetics-solutions-JCrystGrowth3-377.pdf" || true
[ -f "i300005a006.pdf" ] && mv -n -- "i300005a006.pdf" "Oh1982-monolith-converter-transients-IECPDD21-29.pdf" || true
[ -f "1_56.pdf" ] && mv -n -- "1_56.pdf" "Onda1968-gas-liquid-mass-transfer-packed-columns-JCEJ1-56.pdf" || true
[ -f "1-s2.0-S0009250997003850-main.pdf" ] && mv -n -- "1-s2.0-S0009250997003850-main.pdf" "Pan1998-cylindrical-pellet-effectiveness-CES53-933.pdf" || true
[ -f "s11244-018-0948-8.pdf" ] && mv -n -- "s11244-018-0948-8.pdf" "Prins2018-eley-rideal-other-mechanism-TopCatal61-714.pdf" || true
[ -f "Ranz_Marshall_1952_1.pdf" ] && mv -n -- "Ranz_Marshall_1952_1.pdf" "Ranz1952-evaporation-from-drops-pt1-ChemEngProg48-139.pdf" || true
[ -f "Ranz_Marshall_1952_2.pdf" ] && mv -n -- "Ranz_Marshall_1952_2.pdf" "Ranz1952-evaporation-from-drops-pt2-ChemEngProg48-173.pdf" || true
[ -f "1-s2.0-S0263876297800068-main.pdf" ] && mv -n -- "1-s2.0-S0263876297800068-main.pdf" "Richardson1954-sedimentation-fluidisation-pt1-ChERD75-S82-REPRINT1997.pdf" || true
[ -f "1-s2.0-S0376738808003347-main.pdf" ] && mv -n -- "1-s2.0-S0376738808003347-main.pdf" "Robeson2008-upper-bound-revisited-JMembrSci320-390.pdf" || true
[ -f "ie00016a010.pdf" ] && mv -n -- "ie00016a010.pdf" "Rocha1993-structured-packing-hydraulics-pt1-IECR32-641.pdf" || true
[ -f "1-s2.0-S0082078482802816-main.pdf" ] && mv -n -- "1-s2.0-S0082078482802816-main.pdf" "Smith1982-combustion-rates-coal-chars-review-SympCombust19-1045.pdf" || true
[ -f "BF02822675.pdf" ] && mv -n -- "BF02822675.pdf" "Sohn1978-law-of-additive-reaction-times-MetallTransB9B-89.pdf" || true
[ -f "ef00010a006.pdf" ] && mv -n -- "ef00010a006.pdf" "Solomon1988-general-model-coal-devolatilization-FGDVC-EnergyFuels2-405.pdf" || true
[ -f "Syamlal-Rogers-OBrien-1993-MFIX-theory-guide-DOE-METC-94-1004.pdf" ] && mv -n -- "Syamlal-Rogers-OBrien-1993-MFIX-theory-guide-DOE-METC-94-1004.pdf" "Syamlal1993-MFIX-theory-guide-DOE-METC-94-1004.pdf" || true
[ -f "1-s2.0-0009250970850539-main.pdf" ] && mv -n -- "1-s2.0-0009250970850539-main.pdf" "Szekely1970-grain-model-pt1-CES25-1091.pdf" || true
[ -f "rspa.1953.0139.pdf" ] && mv -n -- "rspa.1953.0139.pdf" "Taylor1953-dispersion-solute-in-solvent-ProcRSocA219-186.pdf" || true
[ -f "1-s2.0-S1631074817301091-main.pdf" ] && mv -n -- "1-s2.0-S1631074817301091-main.pdf" "Tayrabekova2018-ethanol-dehydrogenation-copper-CRChimie21-194.pdf" || true
[ -f "Toomey_Johnstone_1952.pdf" ] && mv -n -- "Toomey_Johnstone_1952.pdf" "Toomey1952-gaseous-fluidization-solid-particles-ChemEngProg48-220.pdf" || true
[ -f "1-s2.0-0009250974800898-main.pdf" ] && mv -n -- "1-s2.0-0009250974800898-main.pdf" "Uppal1974-cstr-dynamic-behaviour-CES29-967.pdf" || true
[ -f "Can J Chem Eng - October 1996 - Vanden Bussche - The STAR configuration for methanol synthesis in reversed flow reactors.pdf" ] && mv -n -- "Can J Chem Eng - October 1996 - Vanden Bussche - The STAR configuration for methanol synthesis in reversed flow reactors.pdf" "VandenBussche1996-STAR-reversed-flow-methanol-CJChE74-729.pdf" || true
[ -f "ie50424a010.pdf" ] && mv -n -- "ie50424a010.pdf" "Voorhies1945-carbon-formation-catalytic-cracking-IEC37-318.pdf" || true
[ -f "1-s2.0-0009250962870158-main.pdf" ] && mv -n -- "1-s2.0-0009250962870158-main.pdf" "Wakao1962-random-pore-diffusion-pellets-CES17-825.pdf" || true
[ -f "1-s2.0-0009250978851203-main.pdf" ] && mv -n -- "1-s2.0-0009250978851203-main.pdf" "Wakao1978-particle-to-fluid-transfer-CES33-1375.pdf" || true
[ -f "G1_7.pdf" ] && mv -n -- "G1_7.pdf" "Wammes1991-high-pressure-trickle-bed-hydrodynamics-CET14-406.pdf" || true
[ -f "1-s2.0-0009250956800146-main.pdf" ] && mv -n -- "1-s2.0-0009250956800146-main.pdf" "Wehner1956-boundary-conditions-flow-reactor-CES6-89.pdf" || true
[ -f "1-s2.0-0009250962850052-main.pdf" ] && mv -n -- "1-s2.0-0009250962850052-main.pdf" "Weisz1962-nonisothermal-effectiveness-CES17-265.pdf" || true
[ -f "AIChE Journal - May 1966 - Wen - A generalized method for predicting the minimum fluidization velocity.pdf" ] && mv -n -- "AIChE Journal - May 1966 - Wen - A generalized method for predicting the minimum fluidization velocity.pdf" "Wen1966-minimum-fluidization-velocity-AIChEJ12-610.pdf" || true
[ -f "AIChE Journal - September 1995 - Westerterp - Wave model for longitudinal dispersion  Development of the model (1).pdf" ] && mv -n -- "AIChE Journal - September 1995 - Westerterp - Wave model for longitudinal dispersion  Development of the model (1).pdf" "Westerterp1995-wave-model-longitudinal-dispersion-pt1-AIChEJ41-2013.pdf" || true
[ -f "1-s2.0-0017931062900327-main.pdf" ] && mv -n -- "1-s2.0-0017931062900327-main.pdf" "Whitman1923-two-film-theory-IJHMT5-429-REPRINT1962.pdf" || true
[ -f "1-s2.0-037673889500102I-main.pdf" ] && mv -n -- "1-s2.0-037673889500102I-main.pdf" "Wijmans1995-solution-diffusion-model-JMembrSci107-1.pdf" || true
[ -f "Wilke_1952.pdf" ] && mv -n -- "Wilke_1952.pdf" "Wilke1950-multicomponent-diffusion-mixture-rule-ChemEngProg46-95.pdf" || true
[ -f "AIChE Journal - January 1989 - Xu.pdf" ] && mv -n -- "AIChE Journal - January 1989 - Xu.pdf" "Xu1989-methane-steam-reforming-kinetics-AIChEJ35-88.pdf" || true
[ -f "1-s2.0-S0082078455800331-main.pdf" ] && mv -n -- "1-s2.0-S0082078455800331-main.pdf" "Yagi1955-combustion-carbon-particles-fluidized-beds-SympCombust5-231.pdf" || true
[ -f "AIChE Journal - September 1957 - Yagi - Studies on effective thermal conductivities in packed beds.pdf" ] && mv -n -- "AIChE Journal - September 1957 - Yagi - Studies on effective thermal conductivities in packed beds.pdf" "Yagi1957-effective-thermal-conductivities-packed-beds-AIChEJ3-373.pdf" || true
[ -f "1-s2.0-0009250959800683-main.pdf" ] && mv -n -- "1-s2.0-0009250959800683-main.pdf" "Zwietering1959-degree-of-mixing-CES11-1.pdf" || true
[ -f "1-s2.0-0009250956800031-main.pdf" ] && mv -n -- "1-s2.0-0009250956800031-main.pdf" "vanDeemter1956-plate-height-chromatography-CES5-271.pdf" || true
[ -f "1-s2.0-0009250970850734-main.pdf" ] && mv -n -- "1-s2.0-0009250970850734-main.pdf" "vanWelsenaere1970-parametric-sensitivity-runaway-CES25-1503.pdf" || true
[ -f "i260071a001.pdf" ] && mv -n -- "i260071a001.pdf" "vantRiet1979-kLa-stirred-vessels-review-IECPDD18-357.pdf" || true

echo "renamed: $(ls -1 *.pdf | wc -l) pdfs in $P; $(ls -1 duplicates/*.pdf 2>/dev/null | wc -l) duplicates parked"
