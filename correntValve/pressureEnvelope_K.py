import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from scripts.hydraulic import hydr


# Исходные данные
L = 100 # Длина участка, км 
dx = 1 # Шаг расчета параметров, км
c = 1000 # Скорость распространения волн давления в трубопроводе(для нефти примерное значение, зависит от многих факторов), м/с
d = 1000 # Внутренний диаметр, мм
D = 1020 # Наружный диаметр, мм
rho = 1000 # Плотность жидкости, кг/м^3
Ny = 1 # Кинематическая вязкость, сСт
w0 = 2.7 # Скорость жидкости в трубопроводе в стационарном режиме, м/с
a = 220 # Коэффициент аппроксимации характеристики насоса, м
b = 8.88e-07 # Коэффициент аппроксимации характеристики насоса, ч^2/м^5
hp = 50 # Подпор в начале участка, м
tao = 200 # Период времени для расчета параметров жидкости вы трубопроводе, с
pressureRiseTime = 10 # Время, в которое происходит повышение давления, начиная с начала отсчета, с
deltaP = 2.5 # Величина повышения давления, МПа
xsi = 0.21 # Коэффициент местного сопротивления клапана
valveCoord = 80 # Координата установки предохранительного клапана, км 
PstartOpening = 2.5
PendOpening = 1.2 * PstartOpening

fig = go.Figure()
fig1 = go.Figure()


for dValve in [0.01, 0.1, 0.15, 0.2, 0.3, 0.5]:
    Svalve = np.pi * dValve**2 / 4
    
    hydr1 = hydr(d=d,
            D=D,
            Ny=Ny,
            rho=rho,
            dx=dx,
            c=c,
            a=a,
            b=b,
            hp=hp,
            xsi=xsi,
            Svalve=Svalve,
            PstartOpening=PstartOpening,
            PendOpening=PendOpening
            )
    
    
    # Временной интервал измерения, с
    dt = dx * 1000 / c

    Q = w0 * np.pi * d**2 / 4 / 10**6 * 3600 # Расход в м^3/ч
    P0 = rho * 9.81 * (hp+3*(a-b*Q**2)) / 10**6 # Начальное давление, МПа

    # Задаем таблицы скорости и давления
    flowRate = pd.DataFrame(columns=range(0, int(L*1000) + 1, int(dx*1000)), index=[0], data=w0)

    pressure = pd.DataFrame(columns=range(0, int(L*1000) + 1, int(dx*1000)), index=[0], data=P0)

    # Расчет падения давления по длине трубопровода
    for i, x in enumerate(pressure.columns[1:]):
        pressure.iloc[0, i+1] = pressure.iloc[0, 0] - hydr1.i(w0)*x/10**6
        
    pressure = pressure.round(3)

    # Количество строк для добавления
    strCount = int(np.ceil(tao / dt))

    # Границы участков
    leftBound = 0
    rightBound = len(pressure.columns) - 1
    
    
    traceListPressure = [
        go.Scatter(visible=True, x=pressure.columns/1000, y=pressure.iloc[0, :], mode='lines', name=0)
    ]

    traceListRate = [
        go.Scatter(visible=True, x=flowRate.columns/1000, y=flowRate.iloc[0, :], mode='lines', name=0)
    ]

    KvalveLst = [] # Список степеней открытия клапана в каждый момент времени

    # Расчет давления и скорости для указанного времени с заданным интервалом по длине трубопровода
    for counter in range(1, strCount + 1):
        pressure.loc[counter * dt, :] = 0
        flowRate.loc[counter * dt, :] = 0
        
        for i, x in enumerate(pressure.columns):
            
            # Левая граница
            
            if i == leftBound:
                # Параметры рассматриваемой точки
                pressure.iloc[-1, i] = hydr1.Pc(Pb=pressure.iloc[-2, i + 1], wb=flowRate.iloc[-2, i + 1], leftBound=True)
                flowRate.iloc[-1, i] = hydr1.wc(Pb=pressure.iloc[-2, i + 1], wb=flowRate.iloc[-2, i + 1], leftBound=True)
                
            # Правая граница
            
            elif i == rightBound: 
                # Параметры рассматриваемой точки
                
                # Учитываем время изменения давления
                if pressure.index[-1] == pressureRiseTime:
                    Pc = pressure.iloc[-1, i] = pressure.iloc[-2, i] + deltaP
                else:
                    Pc = pressure.iloc[-1, i] = pressure.iloc[-2, i]
                
                flowRate.iloc[-1, i] = hydr1.wc(Pa=pressure.iloc[-2, i - 1], wa=flowRate.iloc[-2, i - 1], rightBound=True, Pc=Pc)
            
            # Предохранительный клапан
            
            # Предыдущая точка
            
            elif (x == (valveCoord-dx) * 10**3) & (counter >= 2):
                
                # Параметры следующей точки
                Pb = pressure.iloc[-2, i + 1]
                wb = hydr1.wc(pressure.iloc[-3, i],
                            pressure.iloc[-3, i + 2],
                            flowRate.iloc[-3, i],
                            flowRate.iloc[-3, i + 2],
                            valve=True
                            )[0]
                
                # Параметры предыдущей точки
                Pa = pressure.iloc[-2, i - 1]
                wa = flowRate.iloc[-2, i - 1]
                    
                pressure.iloc[-1, i] = hydr1.Pc(Pa, Pb, wa, wb)
                flowRate.iloc[-1, i] = hydr1.wc(Pa, Pb, wa, wb)

            # Следующая точка
            
            elif (x == (valveCoord+dx) * 10**3) & (counter >= 2):
                
                # Параметры предыдущей точки
                Pa = pressure.iloc[-2, i - 1]
                wa = hydr1.wc(pressure.iloc[-3, i - 2],
                            pressure.iloc[-3, i],
                            flowRate.iloc[-3, i - 2],
                            flowRate.iloc[-3, i],
                            valve=True
                            )[1]
                
                # Параметры следующей точки
                Pb = pressure.iloc[-2, i + 1]
                wb = flowRate.iloc[-2, i + 1]
                    
                pressure.iloc[-1, i] = hydr1.Pc(Pa, Pb, wa, wb)
                flowRate.iloc[-1, i] = hydr1.wc(Pa, Pb, wa, wb)
                    
            # Сам клапан
            
            elif x == valveCoord * 10**3:
                # Параметры следующей точки
                Pb = pressure.iloc[-2, i + 1]
                wb = flowRate.iloc[-2, i + 1]
                
                # Параметры предыдущей точки
                Pa = pressure.iloc[-2, i - 1]
                wa = flowRate.iloc[-2, i - 1]

                # Параметры рассматриваемой точки
                Pc = pressure.iloc[-1, i] = hydr1.Pc(Pa, Pb, wa, wb, valve=True)
                wprev, wnext, marker = hydr1.wc(Pa, Pb, wa, wb, valve=True)
                
                flowRate.iloc[-1, i] = wprev # Если клапан закрыт обе переменные одинаковы и равны скорости в этой точке
                
                # Степень открытия клапана
                if Pc - 0.101325 < PstartOpening:
                    KvalveLst.append(0)
                    
                elif Pc - 0.101325 > PendOpening:
                    KvalveLst.append(1)
                    
                else:
                    eta = (Pc - 0.101325 - PstartOpening) \
                        / (PendOpening - PstartOpening)
                    
                    KvalveLst.append(eta)
                
                
                
            # Промежуточные точки
            
            else:
                # Параметры следующей точки
                Pb = pressure.iloc[-2, i + 1]
                wb = flowRate.iloc[-2, i + 1]
                
                # Параметры предыдущей точки
                Pa = pressure.iloc[-2, i - 1]
                wa = flowRate.iloc[-2, i - 1]

                # Параметры рассматриваемой точки
                pressure.iloc[-1, i] = hydr1.Pc(Pa, Pb, wa, wb)
                flowRate.iloc[-1, i] = hydr1.wc(Pa, Pb, wa, wb)
                
        if marker:
            flowRatey = flowRate.loc[counter * dt, :(valveCoord-dx)*1000].to_list() + [wprev, wnext] + flowRate.loc[counter * dt, (valveCoord+dx)*1000:].to_list()
            flowRatex = flowRate.loc[counter * dt, :(valveCoord-dx)*1000].index.to_list() + [valveCoord*1000, valveCoord*1000] + flowRate.loc[counter * dt, (valveCoord+dx)*1000:].index.to_list()
            
            traceListRate.append(go.Scatter(visible=False, x=np.array(flowRatex)/1000, y=flowRatey, mode='lines', name=counter * dt))
        else:
            traceListRate.append(go.Scatter(visible=False, x=flowRate.columns/1000, y=flowRate.loc[counter * dt, :], mode='lines', name=counter * dt))
            
        traceListPressure.append(go.Scatter(visible=False, x=pressure.columns/1000, y=pressure.loc[counter * dt, :], mode='lines', name=counter * dt))
        
        
    fig.add_trace(
        go.Scatter(
            x=pressure.columns/1000,
            y=pressure.max(axis=0),
            name=str(int(dValve*1000)) + 'мм'
        )
    )
    
    fig1.add_trace(
        go.Scatter(
            x=pressure.index,
            y=np.array(KvalveLst) * 100,
            name=str(int(dValve*1000)) + 'мм'
        )
    )
    
    
fig1.update_layout(
    height=600,
    width=850,
    xaxis_title='$$t, c$$',
    yaxis_title=r'$$\eta, \text{%}$$',
    yaxis_range=[0, 110],
    showlegend=False
)

fig1.update_xaxes(title_font_size=30)
fig1.update_yaxes(title_font_size=30)


fig.update_layout(
    xaxis_title=r'$$x, км$$',
    yaxis_title=r'$$P, МПа$$',
    height=600,
    width=850,
    yaxis_range=[0, 6],
    showlegend=True
)

fig.update_xaxes(title_font_size=30)
fig.update_yaxes(title_font_size=30)


fig.show()
fig1.show()