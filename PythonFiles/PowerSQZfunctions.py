##Mauricio Gamonal & Eugenio Bianchi

#Final Version: September 2024

#Reference: https://arxiv.org/pdf/2410.11812

#We report the functions used to plot the power spectra. 

#PowerqBD computes the N2LO corrections for a quasi-Bunch-Davies vacuum around a pivot scale ks
def PowerqBD(k, ks, epsilons, order):
    pi = np.pi
    C=np.euler_gamma+np.log(2)-2
    log_k_ks = np.log(k/ks)
    
    epsilon1cs=epsilons[0]
    epsilon2cs=epsilons[1]
    epsilon1Hs=epsilons[2]
    epsilon2Hs=epsilons[3]
    epsilon1Zs=epsilons[4]
    epsilon2Zs=epsilons[5]
    
    term0 = 1
    term1 = 1 + ((2 + 3 * C) * epsilon1cs - 2 * epsilon1Hs - 2 * C * epsilon1Hs + C * epsilon1Zs + (3 * epsilon1cs - 2 * epsilon1Hs + epsilon1Zs) * log_k_ks)
    
    term2 =  (
        (-8 + 3 * C + (9 * C**2)/2 + (9 * pi**2)/8) * epsilon1cs**2 
        - 3 * epsilon1Hs**2 + 2 * C * epsilon1Hs**2 + 2 * C**2 * epsilon1Hs**2 + (pi**2 * epsilon1Hs**2)/2 
        + 4 * epsilon1Hs * epsilon1Zs - C * epsilon1Hs * epsilon1Zs - 2 * C**2 * epsilon1Hs * epsilon1Zs 
        - (pi**2 * epsilon1Hs * epsilon1Zs)/2 - epsilon1Zs**2 + (C**2 * epsilon1Zs**2)/2 + (pi**2 * epsilon1Zs**2)/8
        - (1/8) * epsilon1cs * (
            4 * (-20 + 10 * C + 12 * C**2 + 3 * pi**2) * epsilon1Hs 
            - 2 * (-24 + 4 * C + 12 * C**2 + 3 * pi**2) * epsilon1Zs 
            + (16 + 16 * C + 12 * C**2 - pi**2) * epsilon2cs)
        + 2 * epsilon1Hs * epsilon2Hs + 2 * C * epsilon1Hs * epsilon2Hs + C**2 * epsilon1Hs * epsilon2Hs 
        - (pi**2 * epsilon1Hs * epsilon2Hs)/12 
        - (1/2) * C**2 * epsilon1Zs * epsilon2Zs + (pi**2 * epsilon1Zs * epsilon2Zs)/24
        + ((3 + 9 * C) * epsilon1cs**2 + (2 + 4 * C) * epsilon1Hs**2 
           + epsilon1cs * (-(5 + 12 * C) * epsilon1Hs + epsilon1Zs + 6 * C * epsilon1Zs - 2 * epsilon2cs - 3 * C * epsilon2cs) 
           - epsilon1Hs * (epsilon1Zs + 4 * C * epsilon1Zs - 2 * (1 + C) * epsilon2Hs) 
           + C * epsilon1Zs * (epsilon1Zs - epsilon2Zs)) * log_k_ks
        + (1/2) * (9 * epsilon1cs**2 + 4 * epsilon1Hs**2 - 4 * epsilon1Hs * epsilon1Zs + epsilon1Zs**2 
                   - 3 * epsilon1cs * (4 * epsilon1Hs - 2 * epsilon1Zs + epsilon2cs) 
                   + 2 * epsilon1Hs * epsilon2Hs - epsilon1Zs * epsilon2Zs) * log_k_ks**2
    )
    
    NLO = term1
    N2LO = term1 + term2
    
    if order == 0:
        return k/k
    elif order == 1:
        return NLO
    elif order == 2:
        return N2LO

#PowSQZ computes the full N2LO power spectrum for a squeezed vacuum (same coefficients as Mathematica notebook).
def PowSQZ(k, ks, epsilons):
    pi = np.pi
    C=np.euler_gamma+np.log(2)-2
    
    # Logarithmic term Log[k/ks]
    log_ks = np.log(k/ks)
    
    epsilon1cs=epsilons[0]
    epsilon2cs=epsilons[1]
    epsilon1Hs=epsilons[2]
    epsilon2Hs=epsilons[3]
    epsilon1Zs=epsilons[4]
    epsilon2Zs=epsilons[5]

    #Sensitive to Vacuum Choice: We use an example
    kc=2.5e-4
    alpha_Re= 1 - (1/2)*(kc/k)**2
    alpha_Im= -(kc/k)
    beta_Re= -((1/2)*(kc/k)**2)*np.cos((2*k)/kc)
    beta_Im= -((1/2)*(kc/k)**2)*np.sin((2*k)/kc)

    #Main code
    term1_inner = (3 + 9*C) * epsilon1cs**2 + (2 + 4*C) * epsilon1Hs**2 + epsilon1Zs - epsilon1Hs * (2 + epsilon1Zs + 4*C*epsilon1Zs) + epsilon1cs * (3 - (5 + 12*C)*epsilon1Hs + epsilon1Zs + 6*C*epsilon1Zs - 2*epsilon2cs - 3*C*epsilon2cs) + 2*(1 + C)*epsilon1Hs*epsilon2Hs + C*epsilon1Zs*(epsilon1Zs - epsilon2Zs)

    term1_alpha_beta = 24 + 3*(-64 + 12*C*(2 + 3*C) + 9*pi**2)*epsilon1cs**2 + 12*(-6 + 4*C*(1 + C) + pi**2)*epsilon1Hs**2 + 3*epsilon1Zs*(8*C + (-8 + 4*C**2 + pi**2)*epsilon1Zs) + 6*epsilon1cs*(8 + (40 - 6*pi**2)*epsilon1Hs + 3*(-8 + pi**2)*epsilon1Zs + 4*C*(3 - 5*epsilon1Hs + epsilon1Zs) + 12*C**2*(-2*epsilon1Hs + epsilon1Zs)) + 3*(-4*(4 + C*(4 + 3*C)) + pi**2)*epsilon1cs*epsilon2cs + 2*epsilon1Hs*(-6*(4 + 4*C**2*epsilon1Zs + (-8 + pi**2)*epsilon1Zs + 2*C*(2 + epsilon1Zs)) + (24 + 12*C*(2 + C) - pi**2)*epsilon2Hs) + (-12*C**2 + pi**2)*epsilon1Zs*epsilon2Zs

    term1_alpha_beta_2 = 24 + 12*(-16 + 6*C + 9*C**2)*epsilon1cs**2 + 24*(-3 + 2*C*(1 + C))*epsilon1Hs**2 - 24*epsilon1Hs*(2 - 4*epsilon1Zs + C*(2 + epsilon1Zs + 2*C*epsilon1Zs)) + 3*epsilon1cs*(8*(2 + 10*epsilon1Hs - 6*epsilon1Zs + C*(3 - (5 + 6*C)*epsilon1Hs + epsilon1Zs + 3*C*epsilon1Zs)) + (-4*(4 + C*(4 + 3*C)) + pi**2)*epsilon2cs) + 2*(24 + 12*C*(2 + C) - pi**2)*epsilon1Hs*epsilon2Hs + epsilon1Zs*(24*C - 24*epsilon1Zs + 12*C**2*(epsilon1Zs - epsilon2Zs) + pi**2*epsilon2Zs)

    term1 = (1/24) * (24 * pi * alpha_Re * beta_Im * term1_inner + alpha_Im**2 * term1_alpha_beta + alpha_Re**2 * term1_alpha_beta + (beta_Im**2 + beta_Re**2) * term1_alpha_beta - 2*alpha_Re*beta_Re * term1_alpha_beta_2 - 2*alpha_Im*(12*pi*beta_Re*term1_inner + beta_Im*term1_alpha_beta_2))

    term2_inner = 9*epsilon1cs**2 + 4*epsilon1Hs**2 - 3*epsilon1cs*(4*epsilon1Hs - 2*epsilon1Zs + epsilon2cs) + 2*epsilon1Hs*(-2*epsilon1Zs + epsilon2Hs) + epsilon1Zs*(epsilon1Zs - epsilon2Zs)

    term2 = (-alpha_Im*(pi*beta_Re*term2_inner + 2*beta_Im*term1_inner) + pi*alpha_Re*beta_Im*term2_inner + alpha_Im**2*term1_inner + alpha_Re**2*term1_inner - 2*alpha_Re*beta_Re*term1_inner + (beta_Im**2 + beta_Re**2)*term1_inner) * np.log(k/ks)

    term3 = 0.5 * ((alpha_Im - beta_Im)**2 + (alpha_Re - beta_Re)**2) * term2_inner * np.log(k/ks)**2

    return term1 + term2 + term3

#Naive squeezing factor Upsilon(k) = |\alpha-\beta|^2
def sqz_factor(k):
    #Sensitive to Vacuum Choice: We use an example
    kc=2.5e-4
    alpha_Re= 1 - (1/2)*(kc/k)**2
    alpha_Im= -(kc/k)
    beta_Re= -((1/2)*(kc/k)**2)*np.cos((2*k)/kc)
    beta_Im= -((1/2)*(kc/k)**2)*np.sin((2*k)/kc)
    
    upsilon = alpha_Re**2 + alpha_Im**2 - 2*(alpha_Im*beta_Im + alpha_Re*beta_Re) + beta_Re**2 + beta_Im**2
    
    return upsilon
   
#Corrected squeezing factor Upsilon(k) = |\alpha- \exp(i\delta) \beta|^2
def sqz_factorDELTA(k,ns):
    #Sensitive to Vacuum Choice: We use an example
    kc=2.5e-4
    kstar=0.05
    alpha_Re= 1 - (1/2)*(kc/k)**2
    alpha_Im= -(kc/k)
    beta_Re= -((1/2)*(kc/k)**2)*np.cos((2*k)/kc)
    beta_Im= -((1/2)*(kc/k)**2)*np.sin((2*k)/kc)
    
    delta = (np.pi/2)*(ns-1)-((np.pi/4)*(ns-1)**2)*np.log(k/kstar)
    
    upsilon = alpha_Re**2 + alpha_Im**2 - 2*(alpha_Im*beta_Im + alpha_Re*beta_Re)*np.cos(delta) + 2*(alpha_Re*beta_Im-alpha_Im*beta_Re)*np.sin(delta) + beta_Re**2 + beta_Im**2
    
    return upsilon
