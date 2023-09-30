# y_p é um vetor da soma das soluções particulares em t0 e dy_p é um ve-
# tor da derivada da soma das soluções particulares

function calculate_integrationConstants(C_bar, F11, y_p, dy_p, U0, V0)
   
       # Calculates C2
       C2 = (C_bar.-(2*F11))\(V0.-dy_p.+((F11.-C_bar)*(y_p.-U0)))
   
       # Calculates C1
       C1 = @. U0-y_p-C2
   
       return C1, C2
   
   end