# Specyfikacja: gra o ruchu ładunków elektrycznych

Status: **zatwierdzona**. Claude Code: przeczytaj całość, a przed pisaniem kodu przedstaw plan Etapu 1 (struktura repozytorium, moduły, kolejność prac, testy) i poczekaj na akceptację.

## 1. Cel

Gra edukacyjna, która buduje intuicję dotyczącą rzeczywistego ruchu ładunków. Gracz ustawia nieruchome ładunki na siatce tak, żeby ładunki próbne doleciały z punktu A do obszaru B, omijając przeszkody w postaci ładunków ustawionych przez poziom. Priorytetem jest **realizm fizyczny**: bez uproszczeń wzorów poza jawnie opisanymi i kontrolowanymi przybliżeniami.

## 2. Fizyka

- Silnik jest **trójwymiarowy od początku**. Tryb 2D to ten sam silnik (patrz punkt 4).
- Źródła pola: ładunki punktowe o skończonym promieniu. Trafienie cząstki w ładunek oznacza jej utratę (bez sztucznego łagodzenia pola).
- Równanie ruchu relatywistyczne: `dp/dt = q(E + v×B)`, `p = γmv`. Integrator Borisa z adaptacyjnym krokiem czasowym (krok maleje przy dużym polu lub gradiencie).
- Pole statycznych ładunków: dokładne prawo Coulomba (suma po parach), liczone w f64.
- Drabina modeli oddziaływania między cząstkami próbnymi, wybierana automatycznie na podstawie `v/c` i pokazywana graczowi:
  1. Coulomb między parami (Etap 2).
  2. Przybliżenie Darwina, czyli oddziaływanie magnetyczne do rzędu `v²/c²` (później).
  3. Pełne równania Maxwella na siatce (FDTD), z opóźnieniem i promieniowaniem (później, opcjonalnie).
- Punkt B to obszar (detektor), bo w polu elektrostatycznym nie ma stabilnej równowagi (twierdzenie Earnshawa).
- Poziom definiuje energię i kierunek startu cząstek w A.
- Jednostki wewnętrzne bezwymiarowe, z jawnie zapisanym przeliczeniem na SI.

## 3. Rozgrywka

- Ładunki gracza i poziomu są na razie nieruchome i stoją w węzłach siatki. Trajektorie cząstek są ciągłe.
- Poziom ma zalecaną rozdzielczość siatki, na której wygenerowano rozwiązanie. Gracz może ją zagęścić o całkowity czynnik (2x, 3x…), co zachowuje wszystkie stare węzły. Dowolna rozdzielczość jest dostępna w trybie niestandardowym.
- Poziom ma limit ładunków gracza (liczba, znak, dopuszczalne wielkości).
- **Fale jako test wytrzymałości konfiguracji:** 1 cząstka, przerwa, 2 cząstki, przerwa, 4, przerwa, i tak dalej aż do strumienia. Pierwsza fala jest prawie idealnie próbna, a kolejne coraz bardziej odczuwają wzajemne oddziaływanie (a w przyszłości reakcję materiałów). Przerwy oddzielają kolejne, coraz trudniejsze etapy i pozwalają materiałom się zrelaksować.
- Wynik (high score) to liczba cząstek, które dotarły do B.
- Szybka pętla informacji zwrotnej: trajektoria pierwszej fali przelicza się na żywo przy każdej zmianie ustawienia, bez przycisku „start”.

## 4. Tryby 2D i 3D

- Jeden rdzeń fizyki, dwie osobne warstwy prezentacji i sterowania. Tryb 2D to osobne doświadczenie: własny styl graficzny, poziomy, ustawienia generatora i prostsze sterowanie.
- **Decyzja:** fizyka jest zawsze prawdziwie trójwymiarowa, także w trybie 2D. Tryb 2D to przekrój świata 3D. Wszystkie ładunki leżą w płaszczyźnie symetrii, a cząstki startują w niej z prędkością równoległą do niej, więc jej nie opuszczają. Pole jest prawdziwym polem 3D (`1/r²`), więc intuicja przenosi się do trybu 3D. Fizyka „płaska” z polem `1/r` jest wykluczona.
- Tryb 3D: edycja warstwami (bieżący przekrój siatki jest aktywny, pozostałe półprzezroczyste).

## 5. Sterowanie (klawiatura)

- Kursor skacze po węzłach siatki. Strzałki przesuwają go w warstwie, PageUp i PageDown zmieniają warstwę (tylko 3D).
- Klawisze: postaw ładunek, usuń, zmień znak, zmień wielkość.
- Zaznaczanie grupy, kopiowanie, wklejanie w miejscu kursora, odbicie lustrzane grupy względem wybranej płaszczyzny.
- Powtórzenia z liczbą (np. „wklej 5 razy co 2 komórki”).
- Wszystkie skróty konfigurowalne.
- Oprócz tego (w 3D) tryb sterowania warstwowego myszką ze zmianą warstw przez scroll.

## 6. Wizualizacja (na żywo, GPU)

- Linie pola startujące z małej sfery wokół każdego ładunku, z liczbą linii proporcjonalną do ładunku.
- Pole wektorowe (strzałki) na wybranym przekroju, mapa potencjału w 2D, powierzchnie stałego potencjału w 3D.
- Pomoce do nauki: zwolnione tempo, wektor siły działającej na cząstkę, pasek energii kinetycznej i potencjalnej w trakcie lotu, wskaźnik aktywnego modelu fizycznego.

## 7. Generator poziomów

1. Losuje ładunki poziomu, punkty A i B oraz parametry startu.
2. Szuka na siatce ustawienia k ładunków gracza, przy którym cząstka trafia w B (symulowane wyżarzanie lub przeszukiwanie wiązkowe). Funkcja celu to najmniejsza odległość trajektorii od B.
3. Odrzuca poziomy trywialne: bez ładunków gracza cząstka nie trafia, a rozwiązania z k−1 ładunkami nie istnieją lub są rzadkie.
4. Wymaga wielu różnych rozwiązań (wielokrotne wyszukiwanie z różnych punktów startowych), żeby poziom był rozwiązywalny także sposobami nieprzewidzianymi przez generator.
5. Sprawdza odporność: przesunięcie jednego ładunku o jedną komórkę nie może zawsze niszczyć rozwiązania.
6. Z punktów 3–5 wylicza miarę trudności.
7. Od Etapu 2 weryfikuje poziom symulacją całej sekwencji fal (drogie, więc dokładnie tylko dla najlepszych kandydatów).
8. Działa offline (natywnie). Poziomy zapisuje jako JSON: ładunki poziomu, A, B, parametry startu, zalecana siatka, limity, rozwiązanie wzorcowe, metryki.

## 8. Materiały (przyszłe etapy, uwzględnić w architekturze teraz)

- Metale: ekwipotencjalne (relaksacja rzędu 1e-19 s, czyli natychmiastowa), z ładunkiem indukowanym.
- Dielektryki: `∇·(ε∇φ) = -ρ`, ładowanie przez zatrzymane cząstki, rozładowanie w przerwach przez upływność.
- Półprzewodniki: model dryfowo-dyfuzyjny (Poisson z równaniami ciągłości dla elektronów i dziur, ruchliwość, dyfuzja, rekombinacja). Najdroższy etap.
- Nadprzewodniki: równania Londonów i efekt Meissnera, sensowne dopiero gdy w grze jest pole magnetyczne.
- Cząstki próbne lecą w próżni lub kanałach. Materiały reagują jako otoczenie (ruch wewnątrz materiału nie jest balistyczny).
- Metody: metoda elementów brzegowych dla metali i dielektryków, solver Poissona na siatce 3D (metoda wielosiatkowa na GPU) dla półprzewodników.

## 9. Architektura i technologie

- Rdzeń fizyki w **Rust**: WebAssembly w grze, natywny program dla generatora. Ten sam kod i f64, więc rozwiązanie z generatora zawsze działa w grze. Deterministyczna symulacja.
- Interfejs `FieldSolver` oddzielający źródło pola od ruchu cząstek (implementacje: Coulomb, później Darwin, Poisson na siatce, elementy brzegowe, FDTD).
- Cząstki jako tablice (osobno pozycje, prędkości), a nie obiekty, żeby fale i strumienie skalowały się naturalnie.
- Obliczenia na GPU przez WebGPU (`wgpu`).
- Front: TypeScript, Vite, Three.js (3D) i osobny renderer 2D (PixiJS lub Canvas). Alternatywa do rozważenia: Bevy.
- Gra działa w przeglądarce, poziomy można udostępniać linkiem.

## 10. Etapy i kryteria akceptacji

1. **Etap 1 – próżnia, jedna cząstka.** Rdzeń 3D (Coulomb, Boris, adaptacyjny krok), tryb 2D jako przekrój, edytor klawiaturą, trajektoria na żywo, linie pola i mapa potencjału, generator z punktami 1–6 z sekcji 7, zapis i odczyt poziomów. Kryteria: testy fizyki przechodzą, trajektoria przelicza się w czasie jednej klatki przy 50 ładunkach, generator produkuje poziomy spełniające punkty 3–5.
2. **Etap 2 – fale.** Oddziaływanie Coulomba między cząstkami, sekwencja fal z przerwami, wynik, generator weryfikujący fale.
3. **Etap 3 – pełny tryb 3D.** Edycja warstwami, wizualizacja 3D.
4. **Etap 4 – metale i dielektryki statyczne.**
5. **Etap 5 – ładowanie materiałów przez cząstki, relaksacja w przerwach.**
6. **Etap 6 – pole magnetyczne, przybliżenie Darwina, nadprzewodniki.**
7. **Etap 7 – półprzewodniki.** Opcjonalnie pełne FDTD.

## 11. Testy fizyki (obowiązkowe od Etapu 1)

- Zachowanie energii w polu statycznym (błąd względny poniżej ustalonego progu).
- Rozpraszanie Rutherforda: kąt odchylenia zgodny z wzorem analitycznym.
- Orbita w polu jednego ładunku przeciwnego znaku (problem Keplera): okres i zachowanie momentu pędu.
- Granica relatywistyczna: prędkość nigdy nie przekracza c.
- Symetria: cząstka startująca w płaszczyźnie symetrii jej nie opuszcza.
- Determinizm: ten sam poziom daje identyczny wynik w WebAssembly i natywnie.

## 12. Otwarte decyzje

- Rust z Three.js czy Bevy.
- Skala fizyczna świata (rozmiar komórki, rodzaj cząstek próbnych, typowe energie).
