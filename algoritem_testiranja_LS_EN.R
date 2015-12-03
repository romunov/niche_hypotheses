# iščemo funkcionalno redundanco v ekosistemu (functional redundancy) -

# 1. najprej naredimo clustering in režemo tako, da dobimo dve skupini
# 2. naredimo nulti model razlike med klastri (guildi) iz vseh vrst tako, da iz naključnih parov izračunamo
#    evklidske razdalje -> histogram
# 3. za vsako združbo (8) izračunamo razliko med centroidoma klastrov (zelen x),
#    ta delta X narišeš (izračunaš p) na histogram. pričakujemo, da bo na zgornjem koncu histograma
#    to ponovimo za vse združbe (8x)
#
# 4. za vsak klastr:
#    preštejemo število vrst v klastru in glede na število (2, 3) naredimo nulti model - izbiramo 2 ali 3 vrste
#    iz vseh vrst glavnega drevesa
# 5. izračunamo minimum spanning tree (MST) iz evklidskih razdalj oz. evklidsko
#    razdaljo za dve vrsti, in to narišemo na histogram oz. izračunamo p (10 grafov)
