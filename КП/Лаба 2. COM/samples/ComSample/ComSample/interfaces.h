#pragma once

#ifndef _INTERFACES
#define _INTERFACES

#include <windows.h>

//IStats используетс€ дл€ получени€ имени автомобил€, а также 
//дл€ отображени€ параметров состо€ни€ реализуемого объекта
DECLARE_INTERFACE_(IStats, IUnknown)
{
	STDMETHOD(DisplayStats)() PURE;
	STDMETHOD(GetPetName)(BSTR * petName) PURE;
};

//IEngine моделирует поведение двигател€. ћожем ускор€ть его, 
//получать значени€ максимальной и текущей скорости
DECLARE_INTERFACE_(IEngine, IUnknown)
{
	STDMETHOD(SpeedUp)() PURE;
	STDMETHOD(GetMaxSpeed)(int* maxSpeed) PURE;
	STDMETHOD(GetCurSpeed)(int* curSpeed) PURE;
};

//ICreateCar используетс€ дл€ присвоени€ имени реализуемому  
//объекту и задани€ максимальной скорости
DECLARE_INTERFACE_(ICreateCar, IUnknown)
{
	STDMETHOD(SetPetName)(BSTR petName) PURE;
	STDMETHOD(SetMaxSpeed)(int maxSp) PURE;
};

#endif // _INTERFACES
