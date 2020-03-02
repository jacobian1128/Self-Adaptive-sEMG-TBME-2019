#include <stdio.h>
#include "WinConsoleCtrl.h"

void SetConsoleBufSize(SHORT x, SHORT y)
{
	COORD bufSz = { x, y };
	SetConsoleScreenBufferSize(GetStdHandle(STD_OUTPUT_HANDLE), bufSz);
}

void SetConsoleWinSize(SHORT x, SHORT y)
{
	//char buf[100];
	//sprintf_s(buf, "mode con:cols=%d lines=%d", x, y);
	//system(buf);

	CONSOLE_SCREEN_BUFFER_INFO info;
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &info);

	if (x > info.dwSize.X && y > info.dwSize.Y)
		SetConsoleBufSize(x, y);
	else if (x > info.dwSize.X)
		SetConsoleBufSize(x, info.dwSize.Y);
	else if (y > info.dwSize.Y)
		SetConsoleBufSize(info.dwSize.X, y);

	SMALL_RECT rect;
	rect.Left = 0;
	rect.Right = x - 1;
	rect.Top = 0;
	rect.Bottom = y - 1;
	SetConsoleWindowInfo(GetStdHandle(STD_OUTPUT_HANDLE), TRUE, &rect);
}

void SetConsoleSize(SHORT bufX, SHORT bufY, SHORT winX, SHORT winY)
{
	if (winX > bufX) bufX = winX;
	if (winY > bufY) bufY = winY;

	COORD bufSz = { bufX, bufY };
	SetConsoleScreenBufferSize(GetStdHandle(STD_OUTPUT_HANDLE), bufSz);

	SMALL_RECT rect;
	rect.Left = 0;
	rect.Right = winX - 1;
	rect.Top = 0;
	rect.Bottom = winY - 1;
	SetConsoleWindowInfo(GetStdHandle(STD_OUTPUT_HANDLE), TRUE, &rect);
}

void SetCCursorPos(int x, int y)
{
	CONSOLE_SCREEN_BUFFER_INFO info;
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &info);
	if (x < 0) x = info.dwCursorPosition.X;
	if (y < 0) y = info.dwCursorPosition.Y;
	COORD pos = { (short)x, (short)y };
	SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), pos);
}

void ShowCCursor(BOOL flag)
{
	CONSOLE_CURSOR_INFO cursorInfo;

	GetConsoleCursorInfo(GetStdHandle(STD_OUTPUT_HANDLE), &cursorInfo);
	cursorInfo.bVisible = flag;
	SetConsoleCursorInfo(GetStdHandle(STD_OUTPUT_HANDLE), &cursorInfo);
}

void SetCTextColor(int tColor)
{
	CONSOLE_SCREEN_BUFFER_INFO info;
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &info);

	if (tColor < 0) tColor = CC_GRAY;
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), (info.wAttributes & 0xF0) | (tColor & 0x0F));
}

void SetCBkgndColor(int bColor)
{
	CONSOLE_SCREEN_BUFFER_INFO info;
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &info);

	if (bColor < 0) bColor = CC_BLACK;
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), ((bColor & 0x0F) << 4) | (info.wAttributes & 0x0F));
}

void SetCColor(int tColor, int bColor)
{
	if (tColor < 0) tColor = CC_GRAY;
	if (bColor < 0) bColor = CC_BLACK;
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), (bColor << 4) | (tColor & 0x0F));
}

COORD GetConsoleBufSize()
{
	CONSOLE_SCREEN_BUFFER_INFO info;
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &info);
	return info.dwSize;
}

COORD GetConsoleWinSize()
{
	CONSOLE_SCREEN_BUFFER_INFO info;
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &info);
	COORD size;
	size.X = abs(info.srWindow.Right - info.srWindow.Left) + 1;
	size.Y = abs(info.srWindow.Bottom - info.srWindow.Top) + 1;
	return size;
}

COORD GetCCursorPos()
{
	CONSOLE_SCREEN_BUFFER_INFO info;
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &info);
	return COORD(info.dwCursorPosition);
}

int GetCTextColor()
{
	CONSOLE_SCREEN_BUFFER_INFO info;
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &info);
	return info.wAttributes & 0x0F;
}

int GetCBkgndColor()
{
	CONSOLE_SCREEN_BUFFER_INFO info;
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &info);
	return (info.wAttributes & 0xF0) >> 4;
}

void GetCColor(int * tColor, int * bColor)
{
	CONSOLE_SCREEN_BUFFER_INFO info;
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &info);
	*tColor = info.wAttributes & 0x0F;
	*bColor = (info.wAttributes & 0xF0) >> 4;
}

