# output test

    Code
      pillar(I(1:3))
    Output
      <pillar>
      <I<int>>
             1
             2
             3
    Code
      pillar(I(list(1, 1:2, 1:3)))
    Output
      <pillar>
      <I<list>>
      <dbl [1]>
      <int [2]>
      <int [3]>

