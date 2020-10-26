// dllmain.h : Declaration of module class.

class CCarComponentModule : public ATL::CAtlDllModuleT<CCarComponentModule>
{
public :
	DECLARE_LIBID(LIBID_CarComponentLib)
	DECLARE_REGISTRY_APPID_RESOURCEID(IDR_CARCOMPONENT, "{09b937df-9d93-4edb-9677-35ac8800d588}")
};

extern class CCarComponentModule _AtlModule;
