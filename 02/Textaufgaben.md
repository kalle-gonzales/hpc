# Nr. 1
Die ungefähre Größe des L1-Cache würde ich nach dem erkennen von "Treppen" in den Daten auf 2048 KB schätzen. Bei Überschreitung der Grenze erhöht sich die Dauer sprunghaft massiv. Das würde ich durch die länger benötigte Zeit beim nutzen des Next-Level Cache erklären. Ebenso erklärt sich dadurch, dass die Übertragungsgeschwindigkeit in data.dat, bei der Stelle den starken Sprung hat.
Die Herstellerangaben des Prozessors bzw. Die Daten, welche ich über das Terminal erhalte sind die folgenden:


LEVEL1_ICACHE_SIZE                 32768 <br>
LEVEL1_DCACHE_SIZE                 32768 <br>
LEVEL2_CACHE_SIZE                  1048576 <br>
LEVEL3_CACHE_SIZE                  11534336 <br>
LEVEL4_CACHE_SIZE                  0


Mit den angaben habe ich nicht gerechnet.
# Nr. 2
Die Unterschiede in der Geschwindigkeit beim wechsel der Schleifenanordnung, erklären sich dadurch, dass in der Matrix benachbarte Daten in column-major-order gespeichert mitunter weit auseinander liegen können und damit auf verschiedenen Cache-lines liegen. Wenn die Schleifenanordnung geändert wird, dann können die Cachelines "besser" liegen.

Die transponierte Matrix hat hinsichtlich der column-major-order eine bessere Ausnutzung der Cachelines, da die Matrix B column-weise multipliziert wird. Die gesamtzahl der Cache-Zugriffe sinkt damit.