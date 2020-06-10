#ifndef TFcoreMFCDlg_H_
#define TFcoreMFCDlg_H_
// TFcoreMFCDlg.h: 헤더 파일
//

#pragma once

#ifndef TFcore_H_
#include "TFcore.cuh"
#endif

#include <vtkAutoInit.h>

#define vtkRenderingCore_AUTOINIT 3(vtkRenderingOpenGL2,vtkInteractionStyle, vtkRenderingFreeType)
#define vtkRenderingContext2D_AUTOINIT 1(vtkRenderingContextOpenGL2)

#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkChartXY.h>
#include <vtkTable.h>
#include <vtkPlot.h>
#include <vtkFloatArray.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkPen.h>

// CTFcoreMFCDlg 대화 상자
class CTFcoreMFCDlg : public CDialogEx
{
// 생성입니다.
public:
	CTFcoreMFCDlg(CWnd* pParent = nullptr);	// 표준 생성자입니다.

// 대화 상자 데이터입니다.
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_TFCOREMFC_DIALOG };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV 지원입니다.


// 구현입니다.
protected:
	HICON m_hIcon;

	// 생성된 메시지 맵 함수
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()

private:
	TFcore TF;
	static UINT ThreadWorker(LPVOID pParam);
	CWinThread* pThreadWorker;
	HANDLE hThreadWorker;

public:
	void LogPrintf(CString input);
	void LogPrintf(const char* fmt, ...);

public:
	afx_msg void OnBnClickedLoad();

public:
	vtkNew<vtkRenderWindow> m_vtkRenderWindow;

	void InitVTKWindow(void* hWnd);
	void ResizeVTKWindow();
	afx_msg void OnBnClickedVtk();
	afx_msg void OnBnClickedProduct();
	CButton mGPU_product;
	afx_msg void OnBnClickedProcess();
	afx_msg void OnBnClickedInit();
};

#endif