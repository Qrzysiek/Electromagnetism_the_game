# Glass screen

The independent flight does not match the game, and the tolerance was not widened.

The glass is the dielectric box in the level file, ε = 4. The bound charge on equal rectangles uses the deflection-plates rectangle integral. The condition is σ = (λ / 2π) times the outward normal field of everything else, with λ = (ε − 1) / (ε + 1). One constant is added per body so the bound charge of the body and its mirror sums to zero. The condition is averaged on four Gauss points of each rectangle. The ion's image is left out.

Inside the glass, on the line under the fixed charge, the field is 0.23 of the point-charge field. On the launch line, outside the glass, it is 1.4 of the point-charge field. The solve residual is 5.5e-11. Halving the step moves an end state by 2.4e-8. Doubling the panel edge moves an end state by 1.9e-3, and the tolerance was set to ten times that, 1.9e-2, before the comparison.

The reference leaves through the right edge at y = 5.041. The game arrives on the detector at y = 8.481. The ordinate differs by 3.4, the abscissa by 4, and the outcomes differ. With no player charge both flights leave through the bottom, and the abscissa differs by 0.93. All of those gaps are above the tolerance. The same player charge with the glass removed hits that charge.

Edge grading of the same rectangles, tried separately, left the reference end at y = 4.98. The rectangle result has converged, and it is not the game's triangle result. No passing document is committed.
