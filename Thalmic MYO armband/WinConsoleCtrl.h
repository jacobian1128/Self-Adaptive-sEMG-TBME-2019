/*!
\file		WinConsoleCtrl.h
\brief		Header for controling the console window
\author		Dong-hyun Lee, dh_wizard@logonu.com
\date		2017. 09. 27
\version	0.0.1
\copyright	(c) LogonU Corp. All rights reserved.
*/

#pragma once
#include <Windows.h>

/*! Font and background color of a console windows.*/
enum eConsoleColor {
	CC_DEFAULT = -1,	//!< text(gray), background(black).
	CC_BLACK = 0,		//!< Set console text or background color to black.
	CC_DARKBLUE,		//!< Set console text or background color to dark blue.
	CC_DARKGREEN,		//!< Set console text or background color to dark green.
	CC_DARKCYAN,		//!< Set console text or background color to dark cyan.
	CC_DARKRED,			//!< Set console text or background color to dark red.
	CC_DARKMAGENTA,		//!< Set console text or background color to dark magenta.
	CC_DARKYELLOW,		//!< Set console text or background color to dark yellow.
	CC_GRAY,			//!< Set console text or background color to gray.
	CC_DARKGRAY,		//!< Set console text or background color to dark gray.
	CC_BLUE,			//!< Set console text or background color to blue.
	CC_GREEN,			//!< Set console text or background color to green.
	CC_CYAN,			//!< Set console text or background color to cyan.
	CC_RED,				//!< Set console text or background color to red.
	CC_MAGENTA,			//!< Set console text or background color to magenta.
	CC_YELLOW,			//!< Set console text or background color to yellow.
	CC_WHITE			//!< Set console text or background color to white.
};

/*!
\brief Sets the screen buffer size of the console window.
\details This function set the screen buffer size of the console window.
\param x a buffer size of the console window's width.
\param y a buffer size of the console window's height.
*/
void SetConsoleBufSize(SHORT x, SHORT y);

/*!
\brief Sets the console window size.
\details This function set the console window size.
\param x a size of the console window's width.
\param y a size of the console window's height.
*/
void SetConsoleWinSize(SHORT x, SHORT y);

/*!
\brief Set the screen buffer and window size.
\details This function set the screen buffer and console window size.\n
\param bufX a buffer size of the console window's width.
\param bufY a buffer size of the console window's height.
\param winX a size of the console window's width.
\param winY a size of the console window's height.
*/
void SetConsoleSize(SHORT bufX, SHORT bufY, SHORT winX, SHORT winY);

/*!
\brief Sets cursor position of the console window.
\details This function set the cursor position of the console window.
\param x a x-axis position of the console window.
\param y a y-axis position of the console window.
*/
void SetCCursorPos(int x, int y);

/*!
\brief Sets the visibility of the cursor.
\details This function set the visibility of the cursor.
\param flag If the cursor is visible, /a flag is TRUE.
*/
void ShowCCursor(BOOL flag);

/*!
\brief Sets text color of the console window.
\details This function sets the text color only, without changing the background color of the console window.
\param tColor a text color, range: 0 ~ 15.
\sa eConsoleColor.
*/
void SetCTextColor(int tColor = CC_DEFAULT);

/*!
\brief Sets background color of the console window.
\details This function sets the background color only, without changing the text color of the console window.
\param bColor a background color, range: 0 ~15
\sa eConsoleColor.
*/
void SetCBkgndColor(int bColor = CC_DEFAULT);

/*!
\brief Sets text and background color of the console window.
\details This fundtion sets the text and background color of the console window.
\param tColor a text color, range: 0 ~ 15.
\param bColor a background color, range: 0 ~ 15.
\sa eConsoleColor.
*/
void SetCColor(int tColor = CC_DEFAULT, int bColor = CC_DEFAULT);

/*!
\brief Gets the screen buffer size of the console window.
\details This function returns the screen buffer size of the console window.
\return the screen buffer size of COOR Structure type.
*/
COORD GetConsoleBufSize();

/*!
\brief Gets the console window size.
\details This function returns the console window size.
\return the console window size of COOR Structure type.
*/
COORD GetConsoleWinSize();

/*!
\brief Gets the current cursor position of the console window.
\details This function returns the current cursor position in the console window.
\return the cursor position of COOR Structure type.
*/
COORD GetCCursorPos();

/*!
\brief Gets the currently set text color of the console window.
\details This function returns the current text color of the console window.
\return text color.
\sa eConsoleColor.
*/
int GetCTextColor();

/*!
\brief Gets the currently set background color of the console window.
\details This function returns the current background color of the console window.
\return background color.
\sa eConsoleColor.
*/
int GetCBkgndColor();

/*!
\brief Get the currently set text and background color of the console window.
\details This function returns the current text and background color of the console window.
\param tColor text color.
\param bColor background color.
\sa eConsoleColor.
*/
void GetCColor(int *tColor, int *bColor);
