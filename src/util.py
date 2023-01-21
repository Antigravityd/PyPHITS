import parsec as pc

mat_identifiers =  {"ac", "actinium", "ag", "silver", "al", "aluminum", "am", "americium", "ar", "argon", "as", "arsenic", "at", "astatine", "au", "gold", "b", "boron", "ba", "barium", "be", "beryllium", "bh", "bohrium", "bi", "bismuth", "bk", "berkelium", "br", "bromine", "c", "carbon", "ca", "calcium", "cd", "cadmium", "ce", "cerium", "cf", "californium", "cl", "chlorine", "cm", "curium", "cn", "copernicium", "co", "cobalt", "cr", "chromium", "cs", "cesium", "cu", "copper", "db", "dubnium", "ds", "darmstadtium", "dy", "dysprosium", "er", "erbium", "es", "einsteinium", "eu", "europium", "f", "fluorine", "fe", "iron", "fl", "flerovium", "fm", "fermium", "fr", "francium", "ga", "gallium", "gd", "gadolinium", "ge", "germanium", "h", "hydrogen", "he", "helium", "hf", "hafnium", "hg", "mercury", "ho", "holmium", "hs", "hassium", "i", "iodine", "in", "indium", "ir", "iridium", "k", "potassium", "kr", "krypton", "la", "lanthanum", "li", "lithium", "lr", "lawrencium", "lu", "lutetium", "lv", "livermorium", "mc", "moscovium", "md", "mendelevium", "mg", "magnesium", "mn", "manganese", "mo", "molybdenum", "mt", "meitnerium", "n", "nitrogen", "na", "sodium", "nb", "niobium", "nd", "neodymium", "ne", "neon", "nh", "nihonium", "ni", "nickel", "no", "nobelium", "np", "neptunium", "o", "oxygen", "og", "oganesson", "os", "osmium", "p", "phosphorus", "pa", "protactinium", "pb", "lead", "pd", "palladium", "pm", "promethium", "po", "polonium", "pr", "praseodymium", "pt", "platinum", "pu", "plutonium", "ra", "radium", "rb", "rubidium", "re", "rhenium", "rf", "rutherfordium", "rg", "roentgenium", "rh", "rhodium", "rn", "radon", "ru", "ruthenium", "s", "sulfur", "sb", "antimony", "sc", "scandium", "se", "selenium", "sg", "seaborgium", "si", "silicon", "sm", "samarium", "sn", "tin", "sr", "strontium", "ta", "tantalum", "tb", "terbium", "tc", "technetium", "te", "tellurium", "th", "thorium", "ti", "titanium", "tl", "thallium", "tm", "thulium", "ts", "tennnessine", "u", "uranium", "v", "vanadium", "w", "tungsten", "xe", "xenon", "y", "yttrium", "yb", "ytterbium", "zn", "zinc", "zr", "zirconium"}
 
def parse_mat_long(s, index):
    
def condense_set(ints): # takes an iterable of numbers and returns PHITS-optimized string,
                        # e.g [1,2,3,5,6,7] -> { 1 - 3 } { 5 - 7 }
    ints = list(set(ints)).sort()

    result = "( "
    glob_start = None
    for i, e in enumerate(ints):
        if i == 0:
            break
        if e - ints[i-1] > 1:
            if glob_start is not None:
                result += f"\{ {glob_start} - {ints[i-1]} \} "
                glob_start = None
            else:
                if i == len(ints) - 1 or ints[i+1] - e > 1:
                    result += f"{e} "
                else:
                    glob_start = e  
    

    result += ')'

    return result
