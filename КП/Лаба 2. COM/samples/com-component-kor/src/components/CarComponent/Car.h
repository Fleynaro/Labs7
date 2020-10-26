// Car.h : Declaration of the CCar

#pragma once
#include "resource.h"       // main symbols

#include "CarComponent_i.h"

#if defined(_WIN32_WCE) && !defined(_CE_DCOM) && !defined(_CE_ALLOW_SINGLE_THREADED_OBJECTS_IN_MTA)
#error "Single-threaded COM objects are not properly supported on Windows CE platform, such as the Windows Mobile platforms that do not include full DCOM support. Define _CE_ALLOW_SINGLE_THREADED_OBJECTS_IN_MTA to force ATL to support creating single-thread COM object's and allow use of it's single-threaded COM object implementations. The threading model in your rgs file was set to 'Free' as that is the only threading model supported in non DCOM Windows CE platforms."
#endif

using namespace ATL;

const int MAX_SPEED = 1000;
const int MAX_NAME_LENGTH = 256;

// CCar
class ATL_NO_VTABLE CCar :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CCar, &CLSID_Car>,
	public IDispatchImpl<ICar, &IID_ICar, &LIBID_CarComponentLib, /*wMajor =*/ 1, /*wMinor =*/ 0>,
	public IDispatchImpl<IEngine, &IID_IEngine, &LIBID_CarComponentLib, 1, 0>,
	public IDispatchImpl<IStats, &IID_IStats, &LIBID_CarComponentLib, 1, 0>,
	public IDispatchImpl<ICreateCar, &IID_ICreateCar, &LIBID_CarComponentLib, 1, 0>
{
public:
	DECLARE_REGISTRY_RESOURCEID(106)

BEGIN_COM_MAP(CCar)
		COM_INTERFACE_ENTRY(ICar)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(IEngine)
		COM_INTERFACE_ENTRY(ICreateCar)
		COM_INTERFACE_ENTRY(IStats)
		END_COM_MAP()

	DECLARE_PROTECT_FINAL_CONSTRUCT()

	HRESULT FinalConstruct();

	void FinalRelease();

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

private:
	BSTR _petName;
	int _maxSpeed = 0;
	int _curSpeed = 0;

public:
};

OBJECT_ENTRY_AUTO(__uuidof(Car), CCar)
