//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   PlotWindow Class
//
//   This class is designed to handle the background operations of
//   running a graphics window.  To take advantage of the primitives,
//   a subclass should be written to handle the interface to the
//   application.
//
//   Created 11/21/95
//
// Modified 4/11/97 - Added X11 functionality
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <math.h>
#include <limits.h>
#include <strings.h>
#include   <string.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "plotwindow.h"
#include <X11/Xlibint.h>
static int pw_colors[NumberOfColors] = {-1};

#define   OneInch   72.0
#define   HardCopyWidth   7.0*OneInch
#define   HardCopyHeight   7.0*OneInch
#define FullScale   10000.0
#define PSFont "Times-Roman"
#define PSFontScale 150
#if defined(__SC__) || defined(__MWERKS__)
#define  DOTSIZE 6
#endif
#define   BLUE   0x01
#define   GREEN   0x02
#define   RED      0x04

namespace levelset {
        
//inline int min(int i, int j) {return i<j ? i : j;}
#define THEDISPLAY itsWindow->dpy
#define THESCREEN  itsWindow->dpy->screens[itsWindow->screen]
//
//   PlotWindow
//
//   Constructor for PlotWindow.  An actual window is optional if only
//   a Postscript file is desired.
//
//   Created 11/21/95
//
    void PlotWindow::Initialize(bool yes_window)
    {
#if defined(__SC__) || defined(__MWERKS__)
        GrafPtr   aPort;   // temporary port pointer
        int      i;         
   
        itsWindow = NULL;
        tempWindow = NULL;
        SetRect(&windowBounds, 5, 40, 390, 425);
        for (i=0; i<8; ++i) {
            blackPat.pat[i] = 0xff;
            grayPat.pat[i] = i%2 ? 0xaa : 0x55;
        };
        PSFile = NULL;
        shade = 1;
        if (yes_window) {
#if defined(__SC__)
            console_options.top = 40;
            console_options.left = 395;
            console_options.nrows = 34;
            console_options.ncols = 39;
            console_options.pause_atexit = 0;
#endif
//      cshow(stdout);
            itsWindow = NewCWindow(0L, &windowBounds, "\pGraphics", true,
                                   noGrowDocProc, (WindowPtr) -1L, true, 0);
            GetPort(&aPort);
            SetPort(itsWindow);
#ifdef __SC__
            TextFont(monaco);
#else
            TextFont(kFontIDMonaco);
#endif
            TextSize(9);
            colorPat = NewPixPat();
            SetPort(aPort);
        };
#endif
#if defined(unix) || defined(__unix)
        itsWindow = NULL;
        tempWindow = NULL;
        StandardFont = NULL;
        PSFile = NULL;
        shade = 1;
        int i, l;

        if (yes_window) {
            int x = 0;
            int y = 0;
            itsWindow = new xpcb;
            itsWindow->xmin = 0;
            itsWindow->xmax = FullScale;
            itsWindow->ymin = 0;
            itsWindow->ymax = FullScale;
            itsWindow->dpy = XOpenDisplay(NULL);
            if (itsWindow->dpy == NULL) {
                std::cerr << "Error: Unable to open display " << XDisplayName(NULL)
                          << std::endl;
                exit(1);
            }
            itsWindow->screen = THEDISPLAY->default_screen;
            itsWindow->width = min(WWidth,THESCREEN.width);
            itsWindow->height = min(WHeight, THESCREEN.height);

            XSetWindowAttributes xswa;
            xswa.background_pixel = THESCREEN.white_pixel;
            xswa.event_mask = ExposureMask;
            xswa.backing_store = Always;

            int depth = THESCREEN.root_depth;
            //Visual* visual = THESCREEN.root_visual;
            Colormap cmap = THESCREEN.cmap;
            i = 5;
            XVisualInfo vinfo_template;
            while (!XMatchVisualInfo(itsWindow->dpy, itsWindow->screen, depth, i--, 
                                     &vinfo_template)) ;

            XColor exact_def;
            if (pw_colors[0] == -1) {
                for (i=0; i<NumberOfColors; ++i) {
                    exact_def.blue = 65535*i/(NumberOfColors-1);
                    exact_def.red = 65535*i/(NumberOfColors-1);
                    exact_def.green = 65535*i/(NumberOfColors-1);
                    exact_def.flags = DoRed | DoBlue | DoGreen;
                    l = XAllocColor(itsWindow->dpy, cmap, &exact_def);
                    pw_colors[i] = exact_def.pixel;
                }
            }
            xswa.colormap = cmap;
   
            itsWindow->win = XCreateWindow(itsWindow->dpy, 
//                                     DefaultRootWindow(itsWindow->dpy),
                                           THESCREEN.root,
                                           x, y, itsWindow->width, itsWindow->height,
                                           0, THESCREEN.root_depth,
                                           InputOutput, CopyFromParent,
                                           CWBackPixel|CWEventMask|
                                           CWBackingStore|CWColormap,
                                           &xswa);
            XSetWindowColormap(itsWindow->dpy, itsWindow->win, cmap);
   
            XSizeHints hints;
            hints.max_width = THESCREEN.width;
            hints.max_height = THESCREEN.height;
            hints.height = itsWindow->height;
            hints.width = itsWindow->width;
            hints.flags = PMaxSize | PSize;
            XSetNormalHints(itsWindow->dpy, itsWindow->win, &hints);
   
            XGCValues gcvals;
            gcvals.foreground = THESCREEN.black_pixel;
            gcvals.function = GXcopy;
            itsWindow->gc = XCreateGC(itsWindow->dpy, itsWindow->win, 
                                      GCForeground|GCFunction, &gcvals);
            XStoreName(itsWindow->dpy, itsWindow->win, "Graphics");
            XMapWindow(itsWindow->dpy, itsWindow->win);
            XSynchronize(itsWindow->dpy, 1);
#ifdef USE_OLD_XLIB_CALLS
            XAnyEvent ev;
            XNextEvent(itsWindow->dpy, (_XEvent*)(&ev));
#else
            XEvent ev;
            XNextEvent(itsWindow->dpy, &ev);
#endif
        }
#endif
    }

//
//   ~PlotWindow
//
//   Denstructor for PlotWindow.  
//
//   Created 11/21/95
//
    PlotWindow::~PlotWindow(void)
    {
        if (PSFile != NULL) 
            EndFile();
        if (pause_at_exit && itsWindow) {
            char c;
            std::cout << "Press return to continue.\n";
            std::cin.get(c);
        }
        if (itsWindow) {
#if defined(__SC__) || defined(__MWERKS__)
            CloseWindow(itsWindow);
            DisposePtr(Ptr(itsWindow));
#endif
#if defined(unix) || defined(__unix)
            XFreeGC(itsWindow->dpy, itsWindow->gc);
            XCloseDisplay(itsWindow->dpy);
            delete itsWindow;
#endif
        }
    }

//
// PlotWindow::ShowPlotWindow
//
// Show or Hide the Plot Window
//
// Created 12/12/95
//
    void PlotWindow::ShowPlotWindow(bool yes)
    {
#if defined(__SC__) || defined(__MWERKS__)
        if (yes && (itsWindow != NULL)) 
            ShowWindow(itsWindow);
        else
            HideWindow(itsWindow);
#endif
#if defined(unix) || defined(__unix)
        if (yes && (itsWindow != NULL)) 
            XUnmapWindow(itsWindow->dpy, itsWindow->win);
        else
            XMapWindow(itsWindow->dpy, itsWindow->win);
#endif
    }

//
// PlotWindow::PlotSegment
//
//   Draw a line segment where the screen is assumed to be a box
//   of 10000x10000 pixels
//
// Created 11/21/95
//
    void PlotWindow::PlotSegment(int x0, int y0, int x1, int y1, 
                                 unsigned short cshade)
    {
#if defined(__SC__) || defined(__MWERKS__)
        GrafPtr   aPort;
        float      x,y;
        Rect      pRect;
#endif
   
        if (itsWindow != NULL) {
#if defined(__SC__) || defined(__MWERKS__)
            GetPort(&aPort);
            SetPort(itsWindow);
            pRect = itsWindow->portRect;
            x = (float(x0))/FullScale*(float(pRect.right)-float(pRect.left));
            y = (float(y0))/FullScale*(float(pRect.bottom)-float(pRect.top));
            MoveTo(pRect.left+Round(x), pRect.bottom-Round(y));
            x = (float(x1))/FullScale*(float(pRect.right)-float(pRect.left));
            y = (float(y1))/FullScale*(float(pRect.bottom)-float(pRect.top));
            LineTo(pRect.left+Round(x), pRect.bottom-Round(y));
            SetPort(aPort);
#endif
#if defined(unix) || defined(__unix)
            int ix0 = Round(float(itsWindow->width*(x0-itsWindow->xmin))
                            /(itsWindow->xmax-itsWindow->xmin));
            int ix1 = Round(float(itsWindow->width*(x1-itsWindow->xmin))
                            /(itsWindow->xmax-itsWindow->xmin));
            int iy0 = Round(float(itsWindow->height*(y0-itsWindow->ymin))
                            /(itsWindow->ymax-itsWindow->ymin));
            int iy1 = Round(float(itsWindow->height*(y1-itsWindow->ymin))
                            /(itsWindow->ymax-itsWindow->ymin));
            int index = (int)(float(cshade)/USHRT_MAX*float(NumberOfColors));
            index = index >= NumberOfColors ? NumberOfColors-1 : index;
            XSetForeground(itsWindow->dpy, itsWindow->gc, (unsigned long)(pw_colors[index]));
            Drawable d = itsWindow->win; 
            XDrawLine(itsWindow->dpy, d, itsWindow->gc, ix0, itsWindow->height - iy0,
                      ix1, itsWindow->height - iy1);
#endif
        };
        if (PSFile != NULL) {
            *PSFile << x0 << " " << y0 << " " << x1 << " " << y1 << " dl\n";
            PSFile->flush();
        };
        return;
    }

//
//   PlotWindow::PlotText
//
//   Plot a string of text with the lower left corner of the text
//   at the specified location.
//
//   Created 11/21/95
//
    void PlotWindow::PlotText(int x0, int y0, const char *text)
    {
#if defined(__SC__) || defined(__MWERKS__)
        GrafPtr   aPort;
        float      x,y;
        Rect      pRect;
#endif
   
        if (itsWindow != NULL) {
#if defined(__SC__) || defined(__MWERKS__)
            GetPort(&aPort);
            SetPort(itsWindow);
            pRect = itsWindow->portRect;
            x = (float(x0))/FullScale*(float(pRect.right)-float(pRect.left));
            y = (float(y0))/FullScale*(float(pRect.bottom)-float(pRect.top));
            MoveTo(pRect.left+Round(x), pRect.bottom-Round(y));
            DrawString(CtoPstr(text));
            PtoCstr((unsigned char*)text);
            SetPort(aPort);
#endif
#if defined(unix) || defined(__unix)
            int ix0 = Round(float(itsWindow->width*(x0-itsWindow->xmin))
                            /(itsWindow->xmax - itsWindow->xmin));
            int iy0 = Round(float(itsWindow->height*(y0-itsWindow->ymin))
                            /(itsWindow->ymax - itsWindow->ymin));
            Drawable d = itsWindow->win;
            XDrawString(itsWindow->dpy, d, itsWindow->gc, ix0, 
                        itsWindow->height - iy0, text, strlen(text));
#endif
        };
        if (PSFile != NULL) {
            *PSFile << x0 << " " << y0 << " moveto\n";
            *PSFile << "(" << text << ") show\n";
            PSFile->flush();
        };
        return;
    }

#define MARGIN 200
    void PlotWindow::PlotText(PlotWindow::WinPos wp, const char* s)
    {
#if defined(unix) || defined(__unix)
        if (itsWindow != NULL) {
            int x0 = 0;
            int y0 = 0;
            int p = 0;
            switch(wp) {
            case LL: x0 = MARGIN; y0 = MARGIN; break; 
            case LR: 
                p = XTextWidth(XQueryFont(itsWindow->dpy,
                                          XGContextFromGC(itsWindow->gc)),
                               s,strlen(s));
                x0 = int(itsWindow->xmax - MARGIN - itsWindow->xmin
                         - p*(itsWindow->xmax - itsWindow->xmin)/itsWindow->width);
                y0 = MARGIN;
                break;
            case UL: x0 = MARGIN; y0 = int(itsWindow->ymax - MARGIN - 10); break; 
            case UR:
                p = XTextWidth(XQueryFont(itsWindow->dpy,
                                          XGContextFromGC(itsWindow->gc)),
                               s,strlen(s));
                x0 = int(itsWindow->xmax - MARGIN - itsWindow->xmin
                         - p*(itsWindow->xmax - itsWindow->xmin)/itsWindow->width);
                y0 = int(itsWindow->ymax - MARGIN - 10);
                break;
            }
            PlotText(x0, y0, s);
        }
        else if (PSFile != NULL) {
            int x0, y0;
            if (wp == LL || wp == LR) 
                y0 = MARGIN;
            else
                y0 = int(FullScale - MARGIN);
            if (wp == LL || wp == UL) 
                x0 = MARGIN;
            else
                x0 = int(FullScale - MARGIN);
            *PSFile << x0 << " " << y0 << " moveto\n";
            *PSFile << "(" << s << ") " 
                    << (wp == LL || wp == UL ? "show" : "rshow") << '\n';
            PSFile->flush();
        }
#endif
    }

      
//
//   PlotWindow::PlotTriangle
//
//   Plot a triangle with the given coordinates and shade
//
//   Created 11/21/95
//
    void PlotWindow::PlotTriangle(int x0, int y0, int x1, int y1, int x2, int y2,
                                  unsigned short cshade)
    {
#if defined(__SC__) || defined(__MWERKS__)
        GrafPtr      aPort;
        float         x,y;
        Rect         pRect;
        PolyHandle   triangle;
        RGBColor      color;
#endif

        if (itsWindow != NULL) {
#if defined(__SC__) || defined(__MWERKS__)
            GetPort(&aPort);
            SetPort(itsWindow);
            pRect = itsWindow->portRect;
            triangle = OpenPoly();
            x = (float(x0))/FullScale*(float(pRect.right)-float(pRect.left));
            y = (float(y0))/FullScale*(float(pRect.bottom)-float(pRect.top));
            MoveTo(pRect.left+Round(x), pRect.bottom-Round(y));
            x = (float(x1))/FullScale*(float(pRect.right)-float(pRect.left));
            y = (float(y1))/FullScale*(float(pRect.bottom)-float(pRect.top));
            LineTo(pRect.left+Round(x), pRect.bottom-Round(y));
            x = (float(x2))/FullScale*(float(pRect.right)-float(pRect.left));
            y = (float(y2))/FullScale*(float(pRect.bottom)-float(pRect.top));
            LineTo(pRect.left+Round(x), pRect.bottom-Round(y));
            x = (float(x0))/FullScale*(float(pRect.right)-float(pRect.left));
            y = (float(y0))/FullScale*(float(pRect.bottom)-float(pRect.top));
            LineTo(pRect.left+Round(x), pRect.bottom-Round(y));
            ClosePoly();
            color.blue = USHRT_MAX-50;
            color.red = cshade;
            color.green = cshade;
            MakeRGBPat(colorPat, &color);
            FillCPoly(triangle, colorPat);
            KillPoly(triangle);
            SetPort(aPort);
#endif
#if defined(unix) || defined(__unix)
            XPoint pts[3];
            pts[0].x = Round(float(itsWindow->width*(x0-itsWindow->xmin))
                             /(itsWindow->xmax - itsWindow->xmin));
            pts[1].x = Round(float(itsWindow->width*(x1-itsWindow->xmin))
                             /(itsWindow->xmax - itsWindow->xmin));
            pts[2].x = Round(float(itsWindow->width*(x2-itsWindow->xmin))
                             /(itsWindow->xmax - itsWindow->xmin));
            pts[0].y = itsWindow->height 
                - Round(float(itsWindow->height*(y0-itsWindow->ymin))
                        /(itsWindow->ymax - itsWindow->ymin));
            pts[1].y = itsWindow->height 
                - Round(float(itsWindow->height*(y1-itsWindow->ymin))
                        /(itsWindow->ymax - itsWindow->ymin));
            pts[2].y = itsWindow->height 
                - Round(float(itsWindow->height*(y2-itsWindow->ymin))
                        /(itsWindow->ymax - itsWindow->ymin));
            int index = (int)(float(cshade)/USHRT_MAX*float(NumberOfColors));
            index = index >= NumberOfColors ? NumberOfColors-1 : index;
            XSetForeground(itsWindow->dpy, itsWindow->gc, (unsigned long)(pw_colors[index]));
            Drawable d = itsWindow->win;
            XFillPolygon(itsWindow->dpy, d, itsWindow->gc, pts, 3, 
                         Convex, CoordModeOrigin);
#endif
        };
        if (PSFile != NULL) {
            *PSFile << float(cshade)/USHRT_MAX << " ";
            *PSFile << x0 << " " << y0 << " " << x1 << " " << y1;
            *PSFile << " " << x2 << " " << y2 << " dt\n";
            PSFile->flush();
        };
        return;
    }

//
//   PlotWindow::PlotPolygon
//
//   Plot a polygon with the given coordinates and shade
//
//   Created 7/19/97
//
    void PlotWindow::PlotPolygon(int num, int* xp, int* yp, unsigned short cshade)
    {
#if defined(__SC__) || defined(__MWERKS__)
        GrafPtr      aPort;
        float        x,y;
        Rect         pRect;
        PolyHandle   polygon;
        RGBColor     color;
#endif
   
        if (itsWindow != NULL) {
#if defined(__SC__) || defined(__MWERKS__)
            GetPort(&aPort);
            SetPort(itsWindow);
            pRect = itsWindow->portRect;
            polygon = OpenPoly();
            x = (float(xp[0]))/FullScale*(float(pRect.right)-float(pRect.left));
            y = (float(yp[0]))/FullScale*(float(pRect.bottom)-float(pRect.top));
            MoveTo(pRect.left+Round(x), pRect.bottom-Round(y));
            for (int i=1; i<num; ++i) {
                x = (float(xp[i]))/FullScale*(float(pRect.right)-float(pRect.left));
                y = (float(yp[i]))/FullScale*(float(pRect.bottom)-float(pRect.top));
                LineTo(pRect.left+Round(x), pRect.bottom-Round(y));
            }
            x = (float(xp[0]))/FullScale*(float(pRect.right)-float(pRect.left));
            y = (float(yp[0]))/FullScale*(float(pRect.bottom)-float(pRect.top));
            LineTo(pRect.left+Round(x), pRect.bottom-Round(y));
            ClosePoly();
            color.blue = USHRT_MAX-50;
            color.red = cshade;
            color.green = cshade;
            MakeRGBPat(colorPat, &color);
            FillCPoly(polygon, colorPat);
            KillPoly(polygon);
            SetPort(aPort);
#endif
#if defined(unix) || defined(__unix)
            XPoint* pts = new XPoint[num];
            for (int i=0; i<num; ++i) {
                pts[i].x = Round(float(itsWindow->width*(xp[i]-itsWindow->xmin))
                                 /(itsWindow->xmax - itsWindow->xmin));
                pts[i].y = itsWindow->height 
                    - Round(float(itsWindow->height*(yp[i]-itsWindow->ymin))
                            /(itsWindow->ymax - itsWindow->ymin));
            }
            Drawable d = itsWindow->win;
            int index = (int)(float(cshade)/USHRT_MAX*float(NumberOfColors));
            if (index >= NumberOfColors) index = NumberOfColors-1;
            XSetForeground(itsWindow->dpy, itsWindow->gc, (unsigned long)(pw_colors[index]));
            XFillPolygon(itsWindow->dpy, d, itsWindow->gc, pts, num, 
                         Convex, CoordModeOrigin);
            delete[] pts;
#endif
        };
        if (PSFile != NULL) {
            *PSFile << float(cshade)/float(USHRT_MAX);
            for (int i=0; i<num; ++i)
                *PSFile << " " << xp[i] << " " << yp[i];
            *PSFile << " " << num-1 << " dp\n";
            PSFile->flush();
        };
        return;
    }


//
//   PlotWindow::PlotDot
//
//   Plot a dot centered at the specified location and color
//
//   Created 11/21/95
//
    void PlotWindow::PlotDot(int x0, int y0, int color)
    {
#if defined(__SC__) || defined(__MWERKS__)
        // Color choices:
        //
        //   black      = 0;
        //   blue      = 1;
        //   green      = 2;
        //   aqua      = 3;
        //   red      = 4;
        //   purple   = 5;
        //   yellow   = 6;
        //   white      = 7;
        //
    
        GrafPtr      aPort;
        float         x,y;
        Rect         pRect;
        Rect         circle;
        RGBColor      clr;
#endif
   
        if (itsWindow != NULL) {
#if defined(__SC__) || defined(__MWERKS__)
            GetPort(&aPort);
            SetPort(itsWindow);
            pRect = itsWindow->portRect;
            x = (float(x0))/FullScale*(float(pRect.right)-float(pRect.left));
            y = (float(y0))/FullScale*(float(pRect.bottom)-float(pRect.top));
            SetRect(&circle, pRect.left+Round(x)-DOTSIZE, pRect.bottom-Round(y)-DOTSIZE,
                    pRect.left+Round(x)+DOTSIZE, pRect.bottom-Round(y)+DOTSIZE);
            clr.blue = color & BLUE ? USHRT_MAX-50 : 0;
            clr.red = color & RED ? USHRT_MAX-50 : 0;
            clr.green = color & GREEN ? USHRT_MAX-50 : 0;
            MakeRGBPat(colorPat, &clr);
            FillCOval(&circle, colorPat);
            SetPort(aPort);
#endif
#if defined(unix) || defined(__unix)
            int ix0 = Round(float(itsWindow->width*(x0-itsWindow->xmin))
                            /(itsWindow->xmax - itsWindow->xmin));
            int iy0 = Round(float(itsWindow->height*(y0-itsWindow->ymin))
                            /(itsWindow->ymax - itsWindow->ymin));
            Drawable d = itsWindow->win;
            XFillArc(itsWindow->dpy, d, itsWindow->gc, ix0-DOTSIZE/2, itsWindow->height - iy0 - DOTSIZE/2, DOTSIZE, DOTSIZE, 0, 360*64);
#endif
        };
        if (PSFile != NULL) {
        };
        return;
    }   

//
//   PlotWindow::PlotOval
//
//   Plot an oval at (x0,y0) with dimensions (w,h)
//
//   Created 11/21/95
//
    void PlotWindow::PlotOval(int x0, int y0, int w, int h)
    {
#if defined(__SC__) || defined(__MWERKS__)
        GrafPtr   aPort;
        float      l, t, r, b;
        Rect      pRect;
        Rect      oval;
        RGBColor   color;
#endif
   
        if (itsWindow != NULL) {
#if defined(__SC__) || defined(__MWERKS__)
            GetPort(&aPort);
            SetPort(itsWindow);
            pRect = itsWindow->portRect;
            l = (float(x0-w))/FullScale*(float(pRect.right)-float(pRect.left));
            b = (float(y0-h))/FullScale*(float(pRect.bottom)-float(pRect.top));
            r = (float(x0+w))/FullScale*(float(pRect.right)-float(pRect.left));
            t = (float(y0+h))/FullScale*(float(pRect.bottom)-float(pRect.top));
            SetRect(&oval, pRect.left+Round(l), pRect.bottom-Round(t),
                    pRect.left+Round(r), pRect.bottom-Round(b));
            FrameOval(&oval);      
            SetPort(aPort);
#endif
#if defined(unix) || defined(__unix)
            int ix0 = Round(float(itsWindow->width*(x0-itsWindow->xmin))
                            /(itsWindow->xmax - itsWindow->xmin));
            int iy0 = Round(float(itsWindow->height*(y0-itsWindow->ymin))
                            /(itsWindow->ymax - itsWindow->ymin));
            int wid = Round(float(itsWindow->width*w)/(itsWindow->xmax-itsWindow->xmin));
            int hgt = Round(float(itsWindow->height*h)/(itsWindow->ymax-itsWindow->ymin));
            Drawable d = itsWindow->win;
            XDrawArc(itsWindow->dpy, d, itsWindow->gc, ix0-wid/2, 
                     itsWindow->height - iy0 - hgt/2, wid, hgt, 0, 360*64);
#endif
        };
        if (PSFile != NULL) {
            *PSFile << x0 << " " << y0 << " " << w/2 << " 0 360 " 
                    << x0+w/2 << " " << y0 << " dc\n";
            PSFile->flush();
        };
        return;
    }

//
//   PlotWindow::PlotFilledOval
//
//   Plot an oval at (x0,y0) with dimensions (w,h)
//
//   Created 11/21/95
//
    void PlotWindow::PlotFilledOval(int x0, int y0, int w, int h)
    {
#if defined(__SC__) || defined(__MWERKS__)
        GrafPtr   aPort;
        float      l, t, r, b;
        Rect      pRect;
        Rect      oval;
        RGBColor   color;
#endif
   
        if (itsWindow != NULL) {
#if defined(unix) || defined(__unix)
            int ix0 = Round(float(itsWindow->width*(x0-itsWindow->xmin))
                            /(itsWindow->xmax - itsWindow->xmin));
            int iy0 = Round(float(itsWindow->height*(y0-itsWindow->ymin))
                            /(itsWindow->ymax - itsWindow->ymin));
            int wid = Round(float(itsWindow->width*w)/(itsWindow->xmax-itsWindow->xmin));
            int hgt = Round(float(itsWindow->height*h)/(itsWindow->ymax-itsWindow->ymin));
            Drawable d = itsWindow->win;
            XFillArc(itsWindow->dpy, d, itsWindow->gc, ix0-wid/2, 
                     itsWindow->height - iy0 - hgt/2, wid, hgt, 0, 360*64);
#endif
        };
        if (PSFile != NULL) {
            *PSFile << x0 << " " << y0 << " " << w/2 << " 0 360 " 
                    << x0+w/2 << " " << y0 << " dc\n";
            PSFile->flush();
        };
        return;
    }

//
//   PlotWindow::Clear
//
//   Erase a rectangle at (x0,y0) with dimensions (w,h)
//
//   Created 8/10/04
//
    void PlotWindow::Clear(int x0, int y0, int w, int h)
    {   
        if (itsWindow != NULL) {
#if defined(unix) || defined(__unix)
            int ix0 = Round(float(itsWindow->width*(x0-itsWindow->xmin))
                            /(itsWindow->xmax - itsWindow->xmin));
            int iy0 = Round(float(itsWindow->height*(y0-itsWindow->ymin))
                            /(itsWindow->ymax - itsWindow->ymin));
            int wid = Round(float(itsWindow->width*w)/(itsWindow->xmax-itsWindow->xmin));
            int hgt = Round(float(itsWindow->height*h)/(itsWindow->ymax-itsWindow->ymin));
            Drawable d = itsWindow->win;
            XClearArea(itsWindow->dpy, d, ix0, 
                       itsWindow->height - iy0 - hgt, wid, hgt, false);
//        XDrawRectangle(itsWindow->dpy, d, itsWindow->gc, ix0, 
//               itsWindow->height - iy0 - hgt, wid, hgt);
#endif
        };
        if (PSFile != NULL) {
        };
        return;
    }

//
//   PlotWindow::SetColor
//
//   Set the graylevel of the color.  
//
//   Created 11/21/95
//
    void PlotWindow::SetColor(int color)
    {
#if defined(__SC__) || defined(__MWERKS__)
        GrafPtr   aPort;
#endif 
   
        if (itsWindow != NULL) {
#if defined(__SC__) || defined(__MWERKS__)
            GetPort(&aPort);
            SetPort(itsWindow);
            if (shade == 1) {
                PenPat(&blackPat);
            } else {
                PenPat(&grayPat);
            };
            SetPort(aPort);
#endif
#if defined(unix) || defined(__unix)
            int index = (int)(float(color)/USHRT_MAX*float(NumberOfColors));
            index = index >= NumberOfColors ? NumberOfColors-1 : index;
            XSetForeground(itsWindow->dpy, itsWindow->gc, (unsigned long)(pw_colors[index]));
#endif
        };
        if (PSFile != NULL) {
            *PSFile << float(color)/USHRT_MAX << " setgray\n";
            PSFile->flush();
        };
        return;
    }

//
//   PlotWindow::SetWidth
//
//   Set the pixel width of the lines
//
//   Created 11/21/95
//
    void PlotWindow::SetWidth(int width)
    {
#if defined(__SC__) || defined(__MWERKS__)
        GrafPtr   aPort;
#endif
   
        if (itsWindow != NULL) {
#if defined(__SC__) || defined(__MWERKS__)
            GetPort(&aPort);
            SetPort(itsWindow);
            PenSize(width, width);
            SetPort(aPort);
#endif
#if defined(unix) || defined(__unix)
            XSetLineAttributes(itsWindow->dpy, itsWindow->gc, width, LineSolid,
                               CapButt, JoinMiter);
#endif
        };
        if (PSFile != NULL) {
            *PSFile << width/(HardCopyWidth/FullScale)/2. << " setlinewidth\n";
            PSFile->flush();
        };
        return;
    }

//
//   PlotWindow::Clear
//
//   Erase the screen
//
//   Created 11/21/95
//
    void PlotWindow::Clear(const bool pause)
    {
#if defined(__SC__) || defined(__MWERKS__)
        GrafPtr   aPort;
#endif
   
        if (itsWindow!= NULL) {
            if (pause) {
                char c;
                std::cout << "Press return to continue.\n";
                std::cin.get(c);
            }
#if defined(__SC__) || defined(__MWERKS__)
            GetPort(&aPort);
            SetPort(itsWindow);
            EraseRect(&(itsWindow->portRect));
            SetPort(aPort);
#endif
#if defined(unix) || defined(__unix)
            XClearWindow(itsWindow->dpy, itsWindow->win);
#endif
        };
        if (PSFile != NULL) {
            EndFile();
            StartNewFile();
        };
        return;
    }

//
//   PlotWindow::StartNewFile
//
//   Open a new postscript file with base name.  Files will be called name##.ps
//   Returns true if the file is successfully opened.
//
//   Created 11/21/95
//
    bool PlotWindow::StartNewFile(const char *name)
    {
        char ans = false;
#ifndef __MWERKS__
        std::ostringstream ostr;
#else
        ostringstream ostr(filename, 50);
#endif
   
        filecnt = 0;
        strcpy(basename, name);
        ostr << name << std::setw(2) << std::setfill('0') << filecnt << ".ps" << std::ends;
        tempWindow = itsWindow;
        itsWindow = NULL;
        PSFile = new std::ofstream(ostr.str().c_str());
        if(PSFile == NULL) {
            std::cout << "iostartfile: could not open output file\n";
        } else {
            *PSFile << "%!PS" << std::endl;
            *PSFile << "%%BoundingBox: 54 180 554 680" << std::endl;
            *PSFile << "/dl {newpath moveto lineto stroke} def" << std::endl;
            *PSFile << "/dt {gsave newpath moveto lineto lineto closepath"
                    << " lshade setgray fill grestore} def" << std::endl;
            *PSFile << "/dp {gsave 3 1 roll newpath moveto {lineto} repeat"
                    << " closepath lshade setgray fill grestore} def" << std::endl;
            *PSFile << "/rshow {dup stringwidth pop 120 exch sub 0"
                    << " rmoveto show} def" << std::endl;
            *PSFile << "/dc {newpath moveto arc stroke} def" << std::endl;
            *PSFile << "/lzero {0} def" << std::endl;
            *PSFile << "/lstep {1} def" << std::endl;
            *PSFile << "/lshade {lstep mul lzero add} def" << std::endl;
            *PSFile << "/" << PSFont << " findfont" << std::endl;
            *PSFile << PSFontScale << " scalefont" << std::endl << "setfont" << std::endl;
            *PSFile << 0.75*OneInch << " " << 2.5*OneInch << " translate" << std::endl;
            *PSFile << HardCopyWidth/FullScale << " " << HardCopyHeight/FullScale
                    << " scale" << std::endl;
            SetWidth(1);
            ans = true;
        };
        itsWindow = tempWindow;
        return (ans);
    }

//
//   PlotWindow::StartNewFile
//
//   Open next file in the sequence.  Close the current file if it is still open.
//
//   Created 11/21/95
//
    bool PlotWindow::StartNewFile(void)
    {
#ifndef __MWERKS__
        std::ostringstream   ostr;
#else
        ostringstream   ostr(filename, 50);
#endif
   
        if (PSFile != NULL) 
            EndFile();
        ++filecnt;
        ostr << basename << std::setw(2) << std::setfill('0') << filecnt << ".ps" << std::ends;
        PSFile = new std::ofstream(ostr.str().c_str());
        if(PSFile == NULL) {
            std::cout << "PlotWindow::StartNewFile: could not open output file" << std::endl;
            return(false);
        } else {
            *PSFile << "%!PS" << std::endl;
            *PSFile << "%%BoundingBox: 54 180 554 680" << std::endl;
            *PSFile << "/dl {newpath moveto lineto stroke} def" << std::endl;
            *PSFile << "/dt {gsave newpath moveto lineto lineto closepath"
                    << " lshade setgray fill grestore} def" << std::endl;
            *PSFile << "/dp {gsave 3 1 roll newpath moveto {lineto} repeat"
                    << " closepath lshade setgray fill grestore} def" << std::endl;
            *PSFile << "/dc {newpath moveto arc stroke} def" << std::endl;
            *PSFile << "/rshow {dup stringwidth pop 120 exch sub 0" 
                    << " rmoveto show} def" << std::endl;
            *PSFile << "/lzero {0} def" << std::endl;
            *PSFile << "/lstep {1} def" << std::endl;
            *PSFile << "/lshade {lstep mul lzero add} def" << std::endl;
            *PSFile << "/" << PSFont << " findfont" << std::endl;
            *PSFile << PSFontScale << " scalefont" << std::endl << "setfont" << std::endl;
            *PSFile << 0.75*OneInch << " " << 2.5*OneInch << " translate" << std::endl;
            *PSFile << HardCopyWidth/FullScale << " " << HardCopyHeight/FullScale << " scale" << std::endl;
            SetWidth(1);
            return(true);
        };
    }

//
//   PlotWindow::EndFile
//
//   Flush buffer and close postscript output file.
//
    void PlotWindow::EndFile(void)
    {
        *PSFile << "showpage" << std::endl;
        PSFile->close();
        delete PSFile;
        PSFile = NULL;
    }

    void PlotWindow::Comment(char* s)
    {
        *PSFile << "%% Begin Comment\n%% " << s
                << "\n%% End Comment\n";
        PSFile->flush();
    }

    void PlotWindow::Comment(std::istream& s)
    {
        char line[78];
        *PSFile << "%% Begin Comment\n";
        do {
            s.getline(line, sizeof(line));
            *PSFile << "%% " << line << '\n';
        } while (s.good() && !s.eof());
        *PSFile << "%% End Comment\n";
        PSFile->flush();
    }


    XImage* PlotWindow::GetImage(void)
    {
        XWindowAttributes xwa;
        //XColor *colors;
        XGetWindowAttributes(itsWindow->dpy, itsWindow->win, &xwa);
        XImage* image = XGetImage(itsWindow->dpy, itsWindow->win, 0, 0, 
                                  itsWindow->width, itsWindow->height, AllPlanes, ZPixmap);
        return image;
    }

    union swapun {
        unsigned long l;
        unsigned short s;
        unsigned char b[4];
    };

// ignore TIFF output
#if 0
    void PlotWindow::SaveToTiff(char* fname)
    {
        XWindowAttributes xwa;
        XGetWindowAttributes(itsWindow->dpy, itsWindow->win, &xwa);
        XImage* image = XGetImage(itsWindow->dpy, itsWindow->win, 0, 0,
                                  itsWindow->width, itsWindow->height, AllPlanes, ZPixmap);
        
        // getxcolors

        XColor* colors = new XColor[NumberOfColors];
        for (int i=0; i<NumberOfColors; i++) {
            colors[i].pixel = pw_colors[i];
            colors[i].pad = 0;
        }

        XQueryColors(itsWindow->dpy, xwa.colormap, colors, NumberOfColors);
        
        // convertimage
        unsigned char* pic = new unsigned char[image->width * image->height];
        unsigned char* pptr = pic;

        /* load up the colormap */
        unsigned char rmap[256], gmap[256], bmap[256];
        for (int i=0; i<NumberOfColors; i++) {
            rmap[i] = colors[i].red   >> 8;
            gmap[i] = colors[i].green >> 8;
            bmap[i] = colors[i].blue  >> 8;
        }
        
        int bits_used = image->bitmap_unit;
        //int bits_per_pixel = image->bits_per_pixel;
        unsigned long pixmask;

        if (image->bits_per_pixel == 32) 
            pixmask = 0xffffffff;
        else 
            pixmask = (((unsigned long) 1) << image->bits_per_pixel) - 1;

        /* if we're on an lsbfirst machine, or the image came from an lsbfirst
           machine, we should flip the bytes around.  NOTE:  if we're on an
           lsbfirst machine *and* the image came from an lsbfirst machine, 
           *don't* flip bytes, as it should work out */

        /* pity we don't have a logical exclusive-or */
        /* determine byte order of the machine we're running on */
        union swapun sw;
        
        sw.l = 0x10000000;
        int isLsbMachine = (sw.b[3]) ? 1 : 0;
        int flipBytes = ( isLsbMachine && image->byte_order != LSBFirst) ||
            (!isLsbMachine && image->byte_order == LSBFirst);

        unsigned char *bptr;
        unsigned short *sptr;
        unsigned long *lptr;
        unsigned char *lineptr;
        unsigned short sval = 0;
        unsigned char tmpbyte;
        unsigned long lval = 0;
        int bit_shift = 0;

        for (int i=0; i<image->height; i++) {
            lineptr = (unsigned char *) image->data + (i * image->bytes_per_line);
            bptr = ((unsigned char *)           lineptr) - 1;
            sptr = ((unsigned short *) lineptr) - 1;
            lptr = ((unsigned long *)  lineptr) - 1;
            bits_used = image->bitmap_unit;

            for (int j=0; j<image->width; j++) {
      
                /* get the next pixel value from the image data */

                if (bits_used == image->bitmap_unit) {  /* time to move on to next b/s/l */
                    switch (image->bitmap_unit) {
                    case 8:  bptr++;  break;
                    case 16: sptr++;  sval = *sptr;
                        if (flipBytes) {   /* swap short */
                            sw.s = sval;
                            tmpbyte = sw.b[0];
                            sw.b[0] = sw.b[1];
                            sw.b[1] = tmpbyte;
                            sval = sw.s;
                        }
                        break;
                    case 32: lptr++;  lval = *lptr;
                        if (flipBytes) {   /* swap long */
                            sw.l = lval;
                            tmpbyte = sw.b[0];
                            sw.b[0] = sw.b[3];
                            sw.b[3] = tmpbyte;
                            tmpbyte = sw.b[1];
                            sw.b[1] = sw.b[2];
                            sw.b[2] = tmpbyte;
                            lval = sw.l;
                        }
                        break;
                    }
                   
                    bits_used = 0;
                    if (image->bitmap_bit_order == MSBFirst) 
                        bit_shift = image->bitmap_unit - image->bits_per_pixel;
                    else 
                        bit_shift = 0;
                }

                unsigned long pixvalue = 0;
                switch (image->bitmap_unit) {
                case 8:  pixvalue = (*bptr >> bit_shift) & pixmask;  break;
                case 16: pixvalue = ( sval >> bit_shift) & pixmask;  break;
                case 32: pixvalue = ( lval >> bit_shift) & pixmask;  break;
                }

                if (image->bitmap_bit_order == MSBFirst) 
                    bit_shift -= image->bits_per_pixel;
                else 
                    bit_shift += image->bits_per_pixel;
                bits_used += image->bits_per_pixel;

      
                /* okay, we've got the next pixel value in 'pixvalue' */
      
                /* use pixel value as an index into colors array */

                //*pptr++ = pixvalue & 0xff;
                int k, found;
                for (k=0, found=0; k<NumberOfColors && !found; ++k) 
                    found = pixvalue == (unsigned long)(pw_colors[k]);
                *pptr++ = k-1;
            }
        }

        // WriteTIFF
        TIFF *tif;

        tif = TIFFOpen(fname, "w");

        if (tif) {

            TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, image->width);
            TIFFSetField(tif, TIFFTAG_IMAGELENGTH, image->height);
            TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
            TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
            TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
            TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
            TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, image->height);

            /* write the image data */

            unsigned char rgb[256];
            unsigned char *tpic = new unsigned char[image->width * image->height];
            unsigned char *tp = tpic;
            TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
            TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
            for (int i=0; i<NumberOfColors; i++) 
                rgb[i] = (rmap[i]*11 + gmap[i]*16 + bmap[i]*5) >> 5;
            int ii;
            for (ii=0, pptr=pic; ii<image->width*image->height; ii++,pptr++) 
                *tp++ = rgb[*pptr];
            TIFFWriteEncodedStrip(tif, 0, tpic, image->width*image->height);
            delete[] tpic;
            TIFFClose(tif);
        }
        
        delete[] colors;
        delete[] pic;
    }


    void PlotWindow::SaveToTiff(const char* bname, const int num, const int digits)
    {
        char fname[256];
        std::ostringstream s;
        
        s << bname << std::setw(digits) << std::setfill('0') << num << ".tif" << std::ends;
        strcpy(fname, s.str().c_str());
        SaveToTiff(fname);
    }
#endif // ignore TIFF output

}
