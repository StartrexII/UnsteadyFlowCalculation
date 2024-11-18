import math

import pandas as pd
import numpy as np


class hydr():
    
    
    def __init__(self, d: float, D:float,
                 Ny: float, rho: float, dx: float,
                 c: float, a: float, b: float,
                 hp: float, xsi: float, Svalve: float,
                 PstartOpening: float
                 ) -> None:
        """initial data for calculation

        Args:
            D (float): Наружный даметр, мм
            d (float): Внутренний диаметр, мм
            Ny (float): Кинематическая вязкость, сСт
            rho (float): Плотность, кг/м^3
            dx (float): Шаг расчета параметров, км
            c (float): Скорость распространения волн давления в трубопроводе, м/c
            a (float): Коэффициент аппроксимации характеристики насоса, м
            b (float): Коэффициент аппроксимации характеристики насоса, ч^2/м^5
            hp (float): Подпор в начале участка, м
            xsi (float): Коэффициент местного сопротивления клапана
            Svalve (float): Площадь клапана, м^2
            PstartOPening (float): Давление начала открытия предохранительного клапана, МПа
            
        """
        
        self.D = D
        self.d = d
        self.Ny = Ny
        self.rho = rho
        self.dx = dx
        self.c = c
        self.a = a
        self.b = b
        self.hp = hp
        self.Patm = 0.101325 # Атмосферное давление, МПа
        self.S = np.pi * (d/1000)**2 / 4 # Площадь трубопровода, м^2
        self.Kvalve = Svalve / xsi**0.5 # Коэффициент пропускной способности клапана, м^2
        self.PstartOpening = PstartOpening
    
    
    def Re(self, w: float) -> float:
        """Число Рейнольдса для трубопровода с заданными параметрами

        Args:
            w (float): Скорость жидкости в данной точке, м/с
        Returns:
            float: Число Рейнольдса
        """
        d = self.d / 1000
        
        Re = w * d / self.Ny * 10**6
        
        return Re
    
    
    def Lambda(self, w: float) -> float:
        """Определение коэффициента гидравлического сопротивления для нефтепровода

        Args:
            w (float): Скорость жидкости в данной точке, м/с
            
        Returns:
            float: коэффициент гидравлического сопротивления
        """
        # Число Рейнольдса нефти в данных условиях
        Re = self.Re(w)
        
        # Вычисление коэффициента гидравлического сопротивления
        if Re <= 0:
            return 0
        
        if Re < 2000:
            Lambda = 64 / Re
        elif 2000 <= Re < 10000:
            Lambda = 0.3164 / Re**0.25
        elif 10000 <= Re :
            Lambda = 0.11 * (68/Re+0.2/self.d)**0.25
        
        return Lambda
    
    
    def i(self, w: float) -> float:
        """Гидравлический уклон(-tgB)

        Args:
            

        Returns:
            float: коэффициент гидравлического уклона, 1/м
        """
        d = self.d / 1000
        rho = self.rho
        
        i = self.Lambda(w) * w * abs(w) * rho / (2*d)
        
        return i
    
    ### Общий случай расчета
    
    def Jb(self, Pb: float, wb: float) -> float:
        """ Отрицательная характеристика

        Args:
            Pb (float): Давление в следующей точке, МПа
            wb (float): Скорость в следующей точке, м/c
            
        Returns:
            float: МПа
        """
        Jb = Pb - self.rho*self.c*wb/10**6 \
            + self.i(wb)*self.dx/10**3
            
        return Jb
        
    
    
    def Ja(self, Pa: float, wa: float):
        """ Положительная характеристика

        Args:
            Pa (float): Давление в предыдущей точке, МПа
            wa (float): Скорость в предыдущей точке, м/c
            
        Returns:
            float: МПа
        """
        Ja = Pa + self.rho*self.c*wa/10**6 \
            - self.i(wa)*self.dx/10**3
            
        return Ja
    
    
    def Pc(self,
           Pa: float|None = None,
           Pb: float|None = None,
           wa: float|None = None,
           wb: float|None = None,
           leftBound: bool = False,
           valve: bool = False
           ) -> float:
        """Давление в рассчитываемой точке

        Args:
            Pa (float | None, optional): Давление в предыдущей точке, МПа. Defaults to None.
            Pb (float | None, optional): Давление в следующей точке, МПа. Defaults to None.
            wa (float | None, optional): Скорость в предыдущей точке, м/c. Defaults to None.
            wb (float | None, optional): Скорость в следующей точке, м/c. Defaults to None.
            leftBound (bool, optional): Расчет для левой границы. Defaults to False.

        Returns:
            float: МПа
        """
        # Левая граница(насосная станция)
        if leftBound:
            w = self.wc(Pb=Pb, wb=wb, leftBound=True)
            A = self.hp + 3*(self.a-self.b*(w*np.pi*self.d**2*900/10**6)**2)
            Pc = self.rho * 9.81 * A
            
            return (Pc/10**6).round(3)
        
        # Предохранительный клапан
        if valve:
            a = self.Pc(Pa, Pb, wa, wb)
            
            if a - self.Patm > self.PstartOpening:
                Kvalve = self.Kvalve
            else:
                Kvalve = 0
            
            e = Kvalve / (2*self.S)
            
            # Расчет давления
            B = 2 * self.c * e * self.rho**0.5
            C = 2 * (a-self.Patm) * 10**6
            
            x = ((B**2+4*C)**0.5-B) / 2
            
            Pc = x**2/(2*10**6) + self.Patm
            
            return Pc.round(3)
        
        
        return ((self.Ja(Pa, wa)+self.Jb(Pb, wb)) * 0.5).round(3)
    
    
    def wc(self,
           Pa: float|None = None,
           Pb: float|None = None,
           wa: float|None = None,
           wb: float|None = None,
           rightBound: bool = False,
           leftBound: bool = False,
           Pc: float|None = None,
           valve: bool = False
           ) -> float | tuple:
        """Скорость в рассчитываемой точке

        Args:
            Pa (float | None, optional): Давление в предыдущей точке, МПа. Defaults to None.
            Pb (float | None, optional): Давление в следующей точке, МПа. Defaults to None.
            wa (float | None, optional): Скорость в предыдущей точке, м/c. Defaults to None.
            wb (float | None, optional): Скорость в следующей точке, м/c. Defaults to None.
            rightBound (bool, optional): Расчет для правой границы_. Defaults to False.
            leftBound (bool, optional): Расчет для левой границы. Defaults to False.
            Pc (float | None, optional): Заданное значение давления в правой границе, МПа. Defaults to None.

        Returns:
            (float | tuple, optional): м/c
        """
        
        
        if rightBound:
            wc = ((self.Ja(Pa, wa)-Pc) \
                / (self.rho*self.c) * 10**6).round(3)
            
            return wc
        
        if leftBound:
            # Коэффициенты квадратного уравнения для нахождения скорости на левой границе
            A = 3 * self.b * np.pi**2 * self.d**4 * 7.9461 * self.rho * 10**-6
            B = self.rho * self.c
            C = self.Jb(Pb, wb)*10**6 - (self.hp+3*self.a)*self.rho*9.81
            
            wc = (((B**2-4*A*C)**0.5-B) / (2*A)).round(3)
            
            return wc
        
        if valve:
            a = self.Pc(Pa, Pb, wa, wb)
            
            if a - self.Patm > self.PstartOpening:
                Kvalve = self.Kvalve
            else:
                Kvalve = 0
            
            e = Kvalve / (2*self.S)
            
            w = self.wc(Pa, Pb, wa, wb)
            
            A = e * (2*(self.Pc(Pa, Pb, wa, wb, valve=True)-self.Patm)/self.rho*10**6)**0.5
            
            wnext = w - A
            wprev = w + A
            
            marker = 0
            if wnext != wprev:
                marker = 1
            
            return (wprev.round(3), wnext.round(3), marker)
            
        
        
        wc = ((self.Ja(Pa, wa)-self.Jb(Pb, wb)) \
            / (2*self.rho*self.c) * 10**6).round(3)
        
        return wc
    
    
        