#include "interfaces.h"
#include "iid.h"
#include <stdio.h>

#define MAX_NAME_LENGTH 64
#define MAX_SPEED 100

class CoCar : public IEngine, public ICreateCar, public IStats
{
	int m_refCount;
	BSTR	m_petName; // ������������� ����� SysAllocString(), 
	// �������� � ����� ����� SysFreeString()
	int		m_maxSpeed; // ������������ ��������
	int		m_currSpeed; // ������� �������� ����r
public:
	// ����������� � ���������� ����r
	CoCar();

	virtual ~CoCar();
	
	// IUnknown
	STDMETHODIMP QueryInterface(REFIID riid, void** pIFace);
	STDMETHODIMP_(DWORD)AddRef();
	STDMETHODIMP_(DWORD)Release();

	// IEngine
	// ���������� IEngine
	STDMETHODIMP SpeedUp();

	STDMETHODIMP GetMaxSpeed(int* maxSpeed);

	STDMETHODIMP GetCurSpeed(int* curSpeed);

	// IStats
	// ���������� � ����r ���������� � ����� ���������
	STDMETHODIMP DisplayStats();

	// ���������� ������� ����� ����������� ������ BSTR 
	STDMETHODIMP GetPetName(BSTR* petName);


	// ICreateCar
	// ���������� ICreateCar
	STDMETHODIMP SetPetName(BSTR petName);

	STDMETHODIMP SetMaxSpeed(int maxSp);

};
