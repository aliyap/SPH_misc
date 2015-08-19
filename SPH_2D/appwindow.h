#ifndef APPWINDOW_H
#define APPWINDOW_H

#include <QMainWindow>
#include "viewer.h"

namespace Ui {
class AppWindow;
}

class AppWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit AppWindow(QWidget *parent = 0);
    ~AppWindow();

private:
    Ui::AppWindow *ui;
    Viewer* m_viewer;
};

#endif // APPWINDOW_H
