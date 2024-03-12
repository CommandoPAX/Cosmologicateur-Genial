import sys, getopt

def generer_monofonic (argv) :

    opts, args = getopt.getopt(argv,"n:l:m:f:k:s:",["ngrid =","Lbox =","wdmass =","fnl =","kmin =","sigma ="])


    ##### Les -- semblent ne pas marcher

    """
    n -> ngrid
    l -> taille de la boite (Mpc / h)
    m -> masse de la matière noire chaude (= 0 => LCDM)
    f -> fnl
    k -> kmin
    s -> sigma
    (f, k, s = 0 => pas de non gaussianités)
    """


    ngrid = 0
    taille = 0
    wdmass = 0  
    fnl = 0
    kmin = 0
    sigma = 0  


    for name, value in opts:
        if name in ["-n", "--ngrid"]:
            ngrid = value
        if name in ["-l", "--Lbox"]:
            taille = value
        if name in ["-m", "--wdmass"]:
            wdmass = value
        if name in ["-f", "--fnl"]:
            fnl = value
        if name in ["-k", "--kmin"]:
            kmin = value
        if name in ["-s", "--sigma"]:
            sigma = value

    print (ngrid, taille, wdmass, fnl, kmin, sigma)
    print(float(wdmass) != 0)

generer_monofonic(sys.argv[1:])