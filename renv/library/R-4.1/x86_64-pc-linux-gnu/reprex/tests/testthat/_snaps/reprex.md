# reprex() works with code that deals with srcrefs

    [1] "``` r"                                                    
    [2] "utils::getParseData(parse(text = 'a'))"                   
    [3] "#>   line1 col1 line2 col2 id parent  token terminal text"
    [4] "#> 1     1    1     1    1  1      3 SYMBOL     TRUE    a"
    [5] "#> 3     1    1     1    1  3      0   expr    FALSE"     
    [6] "```"                                                      

