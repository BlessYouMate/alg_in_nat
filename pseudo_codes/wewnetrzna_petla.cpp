// liczba wymiarów n
// liczba kredytów (np. 10_000)
credits = 10_000

// punkt startowy
x = [10, 10, ..., 10]
Fx = F(x)
credits -= 1

// kodowanie x do binarnej reprezentacji
x_bin = encode(x)

// główna pętla
while credits > 0:
    found_better = false

    // pętla szukania pierwszego lepszego sąsiada
    while not found_better and credits > 0:
        // wygeneruj losowego sąsiada (np. flip jednego bitu)
        x_bin_prime = get_neighbour(x_bin)
        x_prime = decode(x_bin_prime)

        // oblicz F(x')
        Fx_prime = F(x_prime)
        credits -= 1

        // jeśli lepszy -> zaakceptuj
        if Fx_prime < Fx:
            x_bin = x_bin_prime
            x = x_prime
            Fx = Fx_prime
            found_better = true

    // jeśli nie znaleziono lepszego sąsiada -> koniec
    if not found_better:
        break
