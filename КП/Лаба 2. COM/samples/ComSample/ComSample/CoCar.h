#include "interfaces.h"
#include "iid.h"
#include <stdio.h>

#define MAX_NAME_LENGTH 64
#define MAX_SPEED 100

class CoCar : public IEngine, public ICreateCar, public IStats
{
	int m_refCount;
	BSTR	m_petName; // Инициализация через SysAllocString(), 
	// удаление — через вызов SysFreeString()
	int		m_maxSpeed; // Максимальная скорость
	int		m_currSpeed; // Текущая скорость СоСаr
public:
	// Конструктор и деструктор СоСаr
	CoCar();

	virtual ~CoCar();
	
	// IUnknown
	STDMETHODIMP QueryInterface(REFIID riid, void** pIFace);
	STDMETHODIMP_(DWORD)AddRef();
	STDMETHODIMP_(DWORD)Release();

	// IEngine
	// Реализация IEngine
	STDMETHODIMP SpeedUp();

	STDMETHODIMP GetMaxSpeed(int* maxSpeed);

	STDMETHODIMP GetCurSpeed(int* curSpeed);

	// IStats
	// Информация о СоСаr помещается в блоки сообщений
	STDMETHODIMP DisplayStats();

	// Возвращает клиенту копию внутреннего буфера BSTR 
	STDMETHODIMP GetPetName(BSTR* petName);


	// ICreateCar
	// Реализация ICreateCar
	STDMETHODIMP SetPetName(BSTR petName);

	STDMETHODIMP SetMaxSpeed(int maxSp);

};
