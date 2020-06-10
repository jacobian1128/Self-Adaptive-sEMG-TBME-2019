
// TFcoreMFCDlg.cpp: 구현 파일
//

#include "framework.h"
#include "TFcoreMFC.h"
#include "TFcoreMFCDlg.h"
#include "afxdialogex.h"

#include <time.h>


#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// 응용 프로그램 정보에 사용되는 CAboutDlg 대화 상자입니다.
class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// 대화 상자 데이터입니다.
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_ABOUTBOX };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 지원입니다.

// 구현입니다.
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialogEx(IDD_ABOUTBOX)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()


// CTFcoreMFCDlg 대화 상자

CTFcoreMFCDlg::CTFcoreMFCDlg(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_TFCOREMFC_DIALOG, pParent), TF(*this)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CTFcoreMFCDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_PRODUCT, mGPU_product);
}

BEGIN_MESSAGE_MAP(CTFcoreMFCDlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_INIT, &CTFcoreMFCDlg::OnBnClickedInit)
	ON_BN_CLICKED(IDC_LOAD, &CTFcoreMFCDlg::OnBnClickedLoad)
	ON_BN_CLICKED(IDC_VTK, &CTFcoreMFCDlg::OnBnClickedVtk)
	ON_BN_CLICKED(IDC_PRODUCT, &CTFcoreMFCDlg::OnBnClickedProduct)
END_MESSAGE_MAP()


// CTFcoreMFCDlg 메시지 처리기

BOOL CTFcoreMFCDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// 시스템 메뉴에 "정보..." 메뉴 항목을 추가합니다.

	// IDM_ABOUTBOX는 시스템 명령 범위에 있어야 합니다.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != nullptr)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// 이 대화 상자의 아이콘을 설정합니다.  응용 프로그램의 주 창이 대화 상자가 아닐 경우에는
	//  프레임워크가 이 작업을 자동으로 수행합니다.
	SetIcon(m_hIcon, TRUE);			// 큰 아이콘을 설정합니다.
	SetIcon(m_hIcon, FALSE);		// 작은 아이콘을 설정합니다.

	// TODO: 여기에 추가 초기화 작업을 추가합니다.
	InitVTKWindow(GetDlgItem(IDC_FIGURE1)->GetSafeHwnd());
	ResizeVTKWindow();

	//mGPU_product.SetCheck(true);
	//OnBnClickedProduct();

	TF.loadData("data.csv");
	OnBnClickedInit();


	return TRUE;  // 포커스를 컨트롤에 설정하지 않으면 TRUE를 반환합니다.
}

void CTFcoreMFCDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialogEx::OnSysCommand(nID, lParam);
	}
}

// 대화 상자에 최소화 단추를 추가할 경우 아이콘을 그리려면
//  아래 코드가 필요합니다.  문서/뷰 모델을 사용하는 MFC 애플리케이션의 경우에는
//  프레임워크에서 이 작업을 자동으로 수행합니다.

void CTFcoreMFCDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // 그리기를 위한 디바이스 컨텍스트입니다.

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// 클라이언트 사각형에서 아이콘을 가운데에 맞춥니다.
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// 아이콘을 그립니다.
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

// 사용자가 최소화된 창을 끄는 동안에 커서가 표시되도록 시스템에서
//  이 함수를 호출합니다.
HCURSOR CTFcoreMFCDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}



void CTFcoreMFCDlg::OnBnClickedInit()
{
	// TODO: 여기에 컨트롤 알림 처리기 코드를 추가합니다.

	pThreadWorker = AfxBeginThread(ThreadWorker, this);
	hThreadWorker = pThreadWorker->m_hThread;
}


void CTFcoreMFCDlg::OnBnClickedLoad()
{
	// TODO: 여기에 컨트롤 알림 처리기 코드를 추가합니다.
	CString fileName;
	static TCHAR BASED_CODE szFilter[] = _T("data file (*.csv, *.mat) | *.CSV;*.MAT;*.csv;*.map; |");
	CFileDialog dlg(TRUE, NULL, NULL, OFN_HIDEREADONLY, szFilter);

	if (IDOK == dlg.DoModal())	{
		fileName = dlg.GetFileName();
	//	MessageBox(fileName);
	}

	// convert CString to std::string under unicode
	CT2CA pszConvertedAnsiString(fileName);
	string str(pszConvertedAnsiString);
	

	TF.loadData(str);

	LogPrintf("Load data file %s", str.c_str());
	LogPrintf("Data size %d x %d", TF.getNumSamples(), TF.getNumChannels());
}

UINT CTFcoreMFCDlg::ThreadWorker(LPVOID pParam)
{
	CTFcoreMFCDlg* self = (CTFcoreMFCDlg*)pParam;

	clock_t startTime, endTime;
	startTime = clock();
	self->LogPrintf("Begin process...");
	if (self->TF.isGPU_Product()) self->LogPrintf("GPU::Product enabled");
	while (self->TF.getRemaining() > 0) {
		self->TF.getSample();
		self->TF.proceedIteration();
		self->TF.writeResult();

		//self->LogPrintf("Prob %.3e / %.3e", self->TF.getProbMax(), self->TF.getProbThresh());
		if (self->TF.isRegister())
			self->LogPrintf("Registration %d", self->TF.getNumPatterns());
	}
	endTime = clock();

	self->LogPrintf("Process done! %.2f [min]", ((double)(endTime - startTime)) / (60 * CLOCKS_PER_SEC));
	return 0;
}

void CTFcoreMFCDlg::LogPrintf(CString input)
{
	// read current string and add a line
	CString str;
	GetDlgItemText(IDC_LOG, str);
	if (str.GetLength() != 0) str += "\r\n";
	str += input;
	SetDlgItemText(IDC_LOG, str);

	// move cursor to the last line
	CEdit* e = (CEdit*)GetDlgItem(IDC_LOG);
	e->SetFocus();
	e->SetSel(-1);
}

void CTFcoreMFCDlg::LogPrintf(const char* fmt, ...)
{
	// printf style (variable args) overload function
	va_list args;
	char str[1024];

	va_start(args, fmt);
	vsprintf(str, fmt, args);
	va_end(args);

	CString input(str);
	LogPrintf(input);
}


void CTFcoreMFCDlg::InitVTKWindow(void* hWnd)
{
	vtkNew<vtkRenderWindowInteractor> interactor;

	interactor->SetInteractorStyle(vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New());

	vtkNew<vtkRenderer> renderer;
	renderer->SetBackground(0.1, 0.1, 0.1);

	m_vtkRenderWindow->SetParentId(hWnd);
	m_vtkRenderWindow->SetInteractor(interactor);
	m_vtkRenderWindow->AddRenderer(renderer);
	m_vtkRenderWindow->Render();
}

void CTFcoreMFCDlg::ResizeVTKWindow()
{
	CRect rc;
	GetDlgItem(IDC_FIGURE1)->GetClientRect(rc);
	m_vtkRenderWindow->SetSize(rc.Width(), rc.Height());
}


void CTFcoreMFCDlg::OnBnClickedVtk()
{
	// TODO: 여기에 컨트롤 알림 처리기 코드를 추가합니다.
	// Create a table with some points in it
	vtkNew<vtkTable> table;

	vtkNew<vtkFloatArray> arrX;
	arrX->SetName("X Axis");
	table->AddColumn(arrX);

	vtkNew<vtkFloatArray> arrC;
	arrC->SetName("Cosine");
	table->AddColumn(arrC);

	vtkNew<vtkFloatArray> arrS;
	arrS->SetName("Sine");
	table->AddColumn(arrS);

	// Fill in the table with some example values
	int numPoints = 69;
	double inc = 7.5 / ((double)numPoints - 1);
	table->SetNumberOfRows(numPoints);
	for (int i = 0; i < numPoints; ++i)
	{
		table->SetValue(i, 0, i * inc);
		table->SetValue(i, 1, cos(i * inc));
		table->SetValue(i, 2, sin(i * inc));
	}

	// Set up the view
	vtkNew<vtkContextView> view;
	view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);

	// 추가
	view->SetRenderWindow(m_vtkRenderWindow);

	// Add multiple line plots, setting the colors etc
	vtkNew<vtkChartXY> chart;
	view->GetScene()->AddItem(chart);
	vtkPlot* line = chart->AddPlot(vtkChart::LINE);
	line->SetInputData(table, 0, 1);
	line->SetColor(0, 255, 0, 255);
	line->SetWidth(1.0);
	line = chart->AddPlot(vtkChart::LINE);
	line->SetInputData(table, 0, 2);
	line->SetColor(255, 0, 0, 255);
	line->SetWidth(5.0);

	// For dotted line, the line type can be from 2 to 5 for different dash/dot
	// patterns (see enum in vtkPen containing DASH_LINE, value 2):
#ifndef WIN32
	line->GetPen()->SetLineType(vtkPen::DASH_LINE);
#endif
	// (ifdef-ed out on Windows because DASH_LINE does not work on Windows
	//  machines with built-in Intel HD graphics card...)

	//view->GetRenderWindow()->SetMultiSamples(0);

	// Start interactor
	view->GetRenderWindow()->Render();
	view->GetInteractor()->Initialize();
	//view->GetInteractor()->Start();
}


void CTFcoreMFCDlg::OnBnClickedProduct()
{
	// TODO: 여기에 컨트롤 알림 처리기 코드를 추가합니다.
	if (mGPU_product.GetCheck()) TF.enableGPU_Product();
	else						 TF.disableGPU_Product();
}
