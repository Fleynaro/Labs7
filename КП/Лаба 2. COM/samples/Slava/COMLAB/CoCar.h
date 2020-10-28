#include "interfaces.h"
#include "stdio.h"
#include "iid.h"

const int MAX_SPEED = 1000;
const int MAX_NAME_LENGTH = 256;

class CoCar : public IEngine, public ICreateCar, public IStats
{
private:
	ULONG m_refCount;
	BSTR	m_petName; // Инициализация через SysAllocString(), 
	// удаление — через вызов SysFreeString()
	int		m_maxSpeed; // Максимальная скорость
	int		m_currSpeed; // Текущая скорость СоСаr

public:
	CoCar();
	virtual ~CoCar();

	// IUnknown
	STDMETHODIMP QueryInterface(REFIID riid, void** pIFace);
	STDMETHODIMP_(DWORD)AddRef();
	STDMETHODIMP_(DWORD)Release();

	// IEngine
	STDMETHODIMP SpeedUp();
	STDMETHODIMP GetMaxSpeed(int* maxSpeed);
	STDMETHODIMP GetCurSpeed(int* curSpeed);

	// IStats
	STDMETHODIMP DisplayStats();
	STDMETHODIMP GetPetName(BSTR* petName);

	// ICreateCar
	STDMETHODIMP SetPetName(BSTR petName);
	STDMETHODIMP SetMaxSpeed(int maxSp);
};

