# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 09:39:53 2022

@author: wachn
"""
from math import sqrt, sin, cos, atan, degrees, tan, pi, radians, atan2, asin
import numpy as np

class Transformacje:
    def __init__(self, model: str = "wgs84"):  #model - argument domyslny
        #https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flattening = (self.a - self.b) / self.a
        self.ecc2 = (2 * self.flattening - self.flattening ** 2)
        print(model, self.b)



    
    def xyz2plh(self, X, Y, Z, output = 'dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        lat
            [stopnie dziesiętne] - szerokość geodezyjna
        lon
            [stopnie dziesiętne] - długośc geodezyjna.
        h : TYPE
            [metry] - wysokość elipsoidalna
        output [STR] - optional, defoulf 
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """
        r   = sqrt(X**2 + Y**2)           # promień
        lat_prev = atan(Z / (r * (1 - self.ecc2)))    # pierwsze przybliilizenie
        lat = 0
        while abs(lat_prev - lat) > 0.000001/206265:    
            lat_prev = lat
            N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
            h = r / cos(lat_prev) - N
            lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lon = atan(Y/X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2);
        h = r / cos(lat) - N       
        if output == "dec_degree":
            return degrees(lat), degrees(lon), h 
        elif output == "dms":
            lat = self.deg2dms(degrees(lat))
            lon = self.deg2dms(degrees(lon))
            return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - output format not defined")

   
    
    def  flh2XYZ(self,f,l,h):
        f = radians(f)
        l = radians(l)
        N = self.a/(1-self.ecc2*(sin(f))**2)**(0.5)
        X = (N+h)*cos(f)*cos(l) 
        Y = (N+h)*cos(f)*sin(l) 
        Z = (N*(1-self.ecc2)+h)*sin(f) 
        
        return(X, Y, Z)

    
    def fl2xy(self,f,l,L0):
        
        """"
        
        ta funkcja przelicza współrzędne geodezyjne do układu Gaussa_Krugera
        
        Argumenty:
            f - szerokosć geodezyjna            | float 
            l - długosć geodezyjna              | float
            L0 - wartosc lambda zero            |int
            a, e2 - parametry elipsoidy GRS'80
            
        Wyniki w kolejnocsci:
            xgk - współrzędna x w układzie Gaussa-Krugera    |float
            ygk - współrzędna y w ukłądzie Gaussa-Krugera    |float
            
            """
        f = radians(f)
        l = radians(l)
        L0 = radians(L0)
        b2 = (self.a**2)*(1-self.ecc2)
        ep2 = (self.a**2-b2)/b2
        t = tan(f)
        n2 = ep2*(cos(f)**2)
        N = self.a/(1-self.ecc2*(sin(f))**2)**(0.5)
        
        
        A0 = 1-(self.ecc2/4)-(3/64)*(self.ecc2**2)-(5/256)*(self.ecc2**3);
        A2 = (3/8)*(self.ecc2 + (self.ecc2**2)/4 + (15/128)*(self.ecc2**3));
        A4 = (15/256)*(self.ecc2**2 + 3/4*(self.ecc2**3));
        A6 = (35/3072)*self.ecc2**3;
        si = self.a*(A0*f - A2*sin(2*f) + A4*sin(4*f) - A6*sin(6*f));
        dL = l - L0
        
        xgk = si + (dL**2/2)*N*sin(f)*cos(f)*(1 + (dL**2/12)*cos(f)**2*(5 - t**2 + 9*n2 + 4*n2**2) + (dL**4/360)*cos(f)**4*(61 - 58*t**2 + t**4 + 14*n2 - 58*n2*t**2))
        ygk = dL*N*cos(f)*(1 + (dL**2/6)*cos(f)**2*(1 - t**2 + n2) + (dL**4/120)*cos(f)**4*(5 - 18*t**2 + t**4 + 14*n2 - 58*n2*t**2))
        
        return(xgk,ygk)

    def neu(self, X,Y,Z,Xs,Ys,Zs):
    
            '''postaci NEU lub ENU - użytkownik ma wybór
            Funkcja liczy współrzędne wektora NEU i zwraca je w 
    
            Argumenty:
            ----------
            X_sr - współrzędna referencyjna X           | typ: float lub int
            Y_sr - współrzędna referencyjna Y           | typ: float lub int
            Z_sr - współrzędna referencyjna Z           | typ: float lub int
            X - współrzędna X punktu                    | typ: float lub int
            Y - współrzędna Y punktu                    | typ: float lub int
            Z - współrzędna Z punktu                    | typ: float lub int
    
            Wyniki w kolejności:
            -------
            NEU - lista złożona z 3-elementów: N, E, U | typ: lista watrości typu int lub float
            lub
            ENU - lista złożona z 3-elementów: E, N, U | typ: lista watrości typu int lub float
            
            Dodatkowy opis:
            ---------------
            Funkcja jest przeznaczona do liczenia wartości N, E i U dla pojedynczych punktów (nie dla list)
            Funkcja zamieni wartości współrzędnych prostokątnych referencyjnych na liczby, pod warunkiem, że są to cyfry.
    
            '''
            XYZsr = []
            dane = [Xs, Ys, Zs]
            XYZr = np.array([X, Y, Z])
            fi, lam, h = self.xyz2plh(X, Y, Z)
            fi = radians(fi)
            lam = radians(lam)
            
            for indeks, wiersz in enumerate(dane): #wiersz zawiera Xs, Ys, Zs
                Xsr = wiersz - XYZr
                #print(Xsr)
                XYZsr.append(Xsr)

            XYZsr = np.array(XYZsr)
                
            Rneu = np.array([[-sin(fi)*cos(lam), -sin(lam), cos(fi)*cos(lam)],
                             [-sin(fi)*sin(lam), cos(lam), cos(fi)*sin(lam)],
                             [cos(fi), 0, sin(fi)]])
    
            neusr = []
            
            for wiersz in XYZsr:
                neu = Rneu.T @ wiersz
                #print(neu)
                neusr.append(neu)
                
            neusr = np.array(neusr)
            return neusr




    def u2000(self,xgk,ygk,L0):
        
        """" 
        
        ta funkcja przelicza współrzędne z układu Gaussa-Krugera do układu płaskiego PL-2000
        
        Argumenty:
            xgk - współrzędna X w ukłądzie Gaussa-Krugera     | float
            ygk - współrzędna Y w układzie Gaussa-Krugera     | float
            L0 - lambda 0                                     | int
        Wyniki w kolejnosci:
            x - współrzędna X w ukłądzie płaskim PL-2000      | float
            y - współrzędna Y w układzie płaskim PL-2000      | float
            
        """
        m2000 = 0.999923
        
        x2000 = xgk * m2000
        y2000 = ygk * m2000 + (L0/3)* 1000000 + 500000
        
        
        return(x2000,y2000)



    def u92(self, xgk,ygk):
        
        """" 
        
        ta funkcja przelicza współrzędne z układu Gaussa-Krugera do układu płaskiego PL-1992
        
        Argumenty:
            xgk - współrzędna X w ukłądzie Gaussa-Krugera     | float
            ygk - współrzędna Y w układzie Gaussa-Krugera     | float
            L0 - lambda 0                                     | int
        Wyniki w kolejnosci:
            x - współrzędna X w ukłądzie płaskim PL-1992      | float
            y - współrzędna Y w układzie płaskim PL-1992      | float
            
        """
    
        
        m0 = 0.9993 ;
        
        x92 = xgk * m0 - 5300000 
        y92 = ygk* m0 + 500000 
    
        return(x92,y92)


    def dist_xy(self, xA, yA, xB, yB):
        """
        Wyznaczenie azymutu AB i odległości skośniej pomiedzy punktami AB
        INPUT:
            xA : [float] : współrzędna x ounktu A
            yA : [float] : współrzędna y ounktu A
            xB : [float] : współrzędna x ounktu B
            yB : [float] : współrzędna y ounktu B
        OUTPUT:
            (Az_deg, dist_AB) - krotka dwuelementowa, gdzie:
                Az_deg : [float] : azymut AB w stopniach dziesiętnych
                dist_AB: [float] : odległość AB w jednostkach jak podano współrzędne.
        EXAMPLE:    
            INP: xA =-45.00; yA = 23.82; xB = 67.98; yB = 34.12 
            RUN: az, dist = azimuth_dist_xy(xA, yA, x_B, y_b)
            OUT: 5.209060574544288, 113.44853635018832
        """
        # wyznaczenie przyrostów współrzednych
        dX = xB - xA
        dY = yB - yA 

        # wyznaczenie długości odcinka AB
        dist_AB =  sqrt(dX**2 +dY**2) # meter
        return dist_AB

    def azymut_elewacja(self, fi, lam, h, xs, ys, zs):
        fi = radians(fi)
        lam = radians(lam)
        
        dane_sat = [xs, ys, zs]
        
        N = self.a/(sqrt(1 - self.ecc2 * (sin(fi))**2))
        
        X = (N + h) * cos(fi) * cos(lam)
        Y = (N + h) * cos(fi) * sin(lam)
        Z = (N * (1 - self.ecc2) + h) * sin(fi)
        
        XYZr = np.array([X, Y, Z])
        
        XYZsr = []
        
        for indeks, wiersz in enumerate(dane_sat): #wiersz zawiera Xs, Ys, Zs
            Xsr = wiersz - XYZr
            #print(Xsr)
            XYZsr.append(Xsr)
            
        XYZsr = np.array(XYZsr)
            
        Rneu = np.array([[-sin(fi)*cos(lam), -sin(lam), cos(fi)*cos(lam)],
                         [-sin(fi)*sin(lam), cos(lam), cos(fi)*sin(lam)],
                         [cos(fi), 0, sin(fi)]])

        neusr = []
        
        for wiersz in XYZsr:
            neu = Rneu.T @ wiersz
            #print(neu)
            neusr.append(neu)
            
        neusr = np.array(neusr)
            
        el = []
        az = []
        
        for wiersz in neusr:
            n = wiersz[0]
            e = wiersz[1]
            u = wiersz[2]
            azm = atan2(e, n)
            elw = asin((u)/sqrt(n**2 + e**2 + u**2))
            az.append(azm)
            el.append(elw)   
            azwall = np.array(az)
            elall = np.array(el)


        return azwall, elall


if __name__ == '__main__':
    #utworzenie obiektu 
    #dane XYZ geocentryczne 
    geo = Transformacje(model = "grs80")
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    phi, lam, h = geo.xyz2plh(X, Y, Z)
    print(phi, lam, h)
    X,Y,Z = geo.flh2XYZ(phi,lam,h)
    print(X,Y,Z)
    xgk, ygk = geo.fl2xy(phi, lam, 19)
    print(xgk, ygk)
    #n, e, u = geo.neu(X, Y, Z, 3660000, 1400000, 5000000)
    x92, y92 = geo.u92(xgk, ygk)
    print(x92, y92)
    xgk, ygk = geo.fl2xy(phi, lam, 21)
    print(xgk, ygk)
    x2000, y2000 = geo.u2000(xgk,ygk,21)
    print(x2000, y2000)
    dist = geo.dist_xy(X, Y, 3667898, 1479876)
    print(dist)
    azw,el = geo.azymut_elewacja(phi, lam, h, 3665456, 67282929, 6473829)
    print(ee)
    neusr = geo.neu(X, Y, Z, 3660000, 1400000, 5000000)
    print(neusr)

    
   
    
#### wrzucić wszystkie wyniki do tablicy i potem odwoływać się do konkretnych kolumn