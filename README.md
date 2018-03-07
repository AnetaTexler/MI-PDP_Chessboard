# MI-PDP_Chessboard
Parallel and distributed programming

Uloha JES: jezdec na šachovnici

Vstupní data:
k = přirozené číslo, k>5, reprezentující délku strany šachovnice S o velikosti kxk 
q = přirozené číslo q<k2-1 reprezentující počet rozmístěných figurek na šachovnici S 
F[1..q] = pole souřadnic rozmístěných figurek na šachovnici S 
J = souřadnice jezdce na šachovnici S 

Pravidla a cíl hry:
Na počátku je na čtvercové šachovnici S rozmístěno q figurek a 1 jezdec (kůň). Tomuto rozmístění figurek budeme říkat počáteční konfigurace. Jeden tah je skok jezdce podle šachových pravidel. Cílem hry je odstranit jezdcem všechny figurky s minimálním počtem tahů tak, aby na šachovnici zůstal samotný jezdec.

Výstup algoritmu:
Minimální posloupnost tahů jezdce vedoucí do cíle (tj. vedoucí do stavu, kdy jezdec zůstává sám na šachovnici). Výstupem bude seznam souřadnic políček, kam jezdec táhne s označením hvězdičkou těch políček, kde došlo k odstranění figurky.

Definice:
Funkce val(X) ohodnotí libovolné políčko X šachovnice S následujícím způsobem: obsahuje-li políčko X figurku, vrací funkce val(X) hodnotu 1, jinak vrací hodnotu 0. Úspěšnost tahu jezdce na políčko X můžeme definovat jako funkci:

U(X) = 8*val(X)-vzdalenost_z_X_k_nejblizsi_figurce

Sekvenční algoritmus:
Sekvenční algoritmus je typu BB-DFS s hloubkou prohledávaného prostoru omezenou horní mezí (viz dále). Řešení vždy existuje. 
Přípustný koncový stav je situace, kdy na šachovnici zůstane samotný jezdec.
Cena, kterou minimalizujeme, je počet tahů, kterými se lze do koncového stavu z počáteční konfigurace dostat.
Algoritmus končí, když je cena rovna dolní mezi (viz dále) a nebo když se prohledá celý stavový prostor do hloubky dané horní mezí. 
Doporučení pro implementaci: při implementaci sekvenčního algoritmu provádějte následníky X aktuálního políčka jezdce Y v pořadí klesající funkce úspěšnosti U(X) tahu na políčko X. Tím vlastně přednostně skáčete na políčka, kde můžete okamžitě sebrat figurku a je šance sebrat další figurku na co nejméně tahů.
Dolní mez na cenu řešení je q, neboť v jednom tahu může jezdec sebrat maximálně jednu figurku, a proto řešení nemůžeme nalézt rychleji než v q tazích.

Horní mez na cenu řešení je k*k-1, neboť toto je minimální počet tahů, kdy jezdec navštíví všechna políčka šachovnice S. Tato mez je v některých situacích příliš pesimistická, proto bude doporučená hodnota horní meze uvedena přímo v datovém souboru s počáteční konfigurací.

Prořezávání větví stromu. Nemá smysl hledat řešení ve stejné nebo větší hloubce, než je už nalezené minimum (akt_min). Proto musí platit:

akt_hloubka + (q-uz_sebrano_figurek) < akt_min,

jinak je větev ukončena.

Paralelní algoritmus:
Paralelní algoritmus je typu Master-Slave.
