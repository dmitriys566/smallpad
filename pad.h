#ifndef PAD_H
#define PAD_H

#include "./ui_integrate.h"
#include "./ui_plot.h"
#include "./ui_solve.h"
#include "./ui_check.h"

//#include <gnuplot-iostream.h>

#include <QVector>
#include <QMainWindow>
#include <QSystemTrayIcon>
#include <QCloseEvent>
#include <QAction>
//#include <X11/Xlib.h>
#include <string>
#include "tinyexpr/tinyexpr.h"
#ifdef Q_OS_LINUX
#include "nativeeventfilter.h"
#endif

#include <QObject>
#include <QAbstractNativeEventFilter>
#include <QVector>
#ifdef Q_OS_LINUX
#include <QX11Info>
#else
#include <windows.h>
#include <winuser.h>
#endif

QT_BEGIN_NAMESPACE
namespace Ui { class pad; }
QT_END_NAMESPACE

class pad : public QMainWindow
{
    Q_OBJECT

public:
    pad(QWidget *parent = nullptr);
    ~pad();
private slots:
    /* The slot that will accept the signal from the event
     * Click on the application icon in the system tray
     */
    void iconActivated(QSystemTrayIcon::ActivationReason reason);
public slots:
    void pad_kill();
    void pad_exit();
    void pad_show();
    void pad_clear();
    void pad_about();
    void pad_integrate();
    void pad_eval_integrate();
    void pad_about_exit();
    void pad_key_exit();
    void pad_key_enter();
    void pad_integrate_exit();
    void pad_plot_do();
    void pad_plot_exit();
    void pad_plot_export();
    void clip_scan();
    void pad_load();
    void pad_save();
    void pad_copy();
    void pad_paste();
    void pad_copylast();
    void pad_tray();
    void pad_plot();
    void pad_solve();
    void pad_solve_do();
    void pad_solve_exit();
    void pad_minimize();
    void pad_tutorial();
    void pad_key();
protected:
    void closeEvent(QCloseEvent * event);
public:
    void evaluate();
    bool pad_wait();
    void pad_check();
//    double eval_func(double *x,void *sd);
    QString selectedText();
    void setSelection(QString text);

protected:
    //    void changeEvent(QEvent *e);
    bool eventFilter(QObject* obj, QEvent *event);
#ifdef Q_OS_WINDOWS
    bool nativeEvent(const QByteArray &eventType, void *message, long *result);
#endif

private slots:
    void on_lineEdit_returnPressed();

private:
    QSystemTrayIcon         * trayIcon;
    QClipboard * m_clip;
#ifdef Q_OS_LINUX
    NativeEventFilter *nativeEventFilter;
#endif
    Ui::pad *ui;
    QDialog *m_about;
    QTimer * m_t;
    QDialog *m_key;
    QDialog *m_integrate;
    QDialog *m_plot;
    QDialog *m_solve;
    Ui::integrate *m_ui_int;
    Ui::plot *m_ui_plot;
    Ui::solve *m_ui_solve;
    Ui::check *m_ui_check;
//    QSystemTrayIcon * trayIcon;
    std::string m_expression;
    std::vector<std::string> m_vars_string;
    std::vector<double> m_vars_double;
    std::vector<std::string> m_history;
    std::vector<te_variable> m_vars;
    std::vector<double> m_results;
    te_parser parser;
    std::string m_ex;
    double m_ans;//Результат последнего вычисления
    int m_ptr=-1;
    int count=0;
    bool m_dirty = false;
};
#endif // PAD_H
