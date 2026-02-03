#include "pad.h"
#include "./ui_pad.h"
#include "./ui_about.h"
#include "./ui_plot.h"
#include "./ui_integrate.h"
#include "./ui_check.h"
#include "libdcl.h"
#include <functional>
#include "tinyexpr/tinyexpr.h"
#include <locale>
#include <string>
#include <chrono>
#include <thread>
#include <clocale>
#include <QTextBrowser>
#include <QTextStream>
#ifdef Q_OS_LINUX
#include "nativeeventfilter.h"
#include <unistd.h>
#include <functional>
#include <chrono>
#include <future>
#include <cstdio>
#endif
//#include <libclipboard/include/libclipboard.h>

//#include "tinyexpr/tinyexpr.c"

//#include <X11/Xlib.h>
//#include <X11/keysym.h>
//#include <stdbool.h>
//#include <stdio.h>

#include <QObject>
#include <QMessageBox>
#include <QClipboard>
#include <QAbstractNativeEventFilter>
#include <QVector>
#include <QMimeData>
#include <QThread>
#include <QTimer>
#ifdef Q_OS_LINUX
#include <QX11Info>
#endif
#include <QFileDialog>

typedef struct
{
    char *expression;
}func_params;

void term()
{
    exit(0);
}

uint32_t crc32(unsigned char *buf, size_t len)
{
    uint_least32_t crc_table[256];
    uint_least32_t crc; int i, j;

    for (i = 0; i < 256; i++)
    {
        crc = i;
        for (j = 0; j < 8; j++)
            crc = crc & 1 ? (crc >> 1) ^ 0xEDB88320UL : crc >> 1;

        crc_table[i] = crc;
    };

    crc = 0xFFFFFFFFUL;

    while (len--)
        crc = crc_table[(crc ^ *buf++) & 0xFF] ^ (crc >> 8);

    return crc ^ 0xFFFFFFFFUL;
}


pad::pad(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::pad)
{
    ui->setupUi(this);
    m_clip = QGuiApplication::clipboard();
    connect(ui->actionExit, SIGNAL(triggered(bool)), this, SLOT(pad_exit()));
    connect(ui->actionClear, SIGNAL(triggered(bool)), this, SLOT(pad_clear()));
    connect(ui->actionAbout, SIGNAL(triggered(bool)), this, SLOT(pad_about()));
    connect(ui->actionIntegrate_function, SIGNAL(triggered(bool)), this, SLOT(pad_integrate()));
    connect(ui->actionLoad, SIGNAL(triggered(bool)), this, SLOT(pad_load()));
    connect(ui->actionSave, SIGNAL(triggered(bool)), this, SLOT(pad_save()));
    connect(ui->actionCut, SIGNAL(triggered(bool)), this, SLOT(pad_copy()));
    connect(ui->actionPaste, SIGNAL(triggered(bool)), this, SLOT(pad_paste()));
    connect(ui->actionCopy_last, SIGNAL(triggered(bool)), this, SLOT(pad_copylast()));
    connect(ui->actionPlot_graph, SIGNAL(triggered(bool)), this, SLOT(pad_plot()));
    connect(ui->actionSolve_nonlinear_equation, SIGNAL(triggered(bool)), this, SLOT(pad_solve()));
    connect(ui->actionMinimize_to_tray, SIGNAL(triggered(bool)), this, SLOT(pad_minimize()));
    connect(ui->actionTutorial, SIGNAL(triggered(bool)), this, SLOT(pad_tutorial()));
    connect(ui->actionEnter_license_key, SIGNAL(triggered(bool)), this, SLOT(pad_key()));
    m_t = new QTimer();
    connect(m_t, SIGNAL(timeout()), this, SLOT(pad_kill()));
    pad_check();
    ui->lineEdit->installEventFilter(this);
#ifdef Q_OS_LINUX
    setlocale(LC_ALL, "en_US.UTF-8");
    std::locale::global(std::locale("en_US.UTF-8"));
    //nativeEventFilter = new NativeEventFilter(this);    // Инициализируем фильтр
    //qApp->installNativeEventFilter(nativeEventFilter);
    //connect(nativeEventFilter, &NativeEventFilter::activated, this, &pad::clip_scan);
    //nativeEventFilter->setShortcut();   // Устанавилваем хоткей
#else
    RegisterHotKey((HWND)pad::winId(),   // Set the system identifier of the widget window that will handle the HotKey
                   100,                         // Set identifier HotKey
                   MOD_CONTROL | MOD_ALT,         // Set modifiers
                   VK_OEM_PLUS);                        // We define hotkeys for HotKey
#endif
    trayIcon = new QSystemTrayIcon(this);
    trayIcon->setIcon(QIcon(":/icon_calculator.png"));
    trayIcon->setToolTip("SmallPAD" "\n"
                         "Калькулятор выражений");
    /* After that create a context menu of two items */
    QMenu * menu = new QMenu(this);
    QAction * viewWindow = new QAction(trUtf8("Развернуть окно"), this);
    QAction * quitAction = new QAction(trUtf8("Выход"), this);

    /* connect the signals clicks on menu items to by appropriate slots.
     * The first menu item expands the application from the tray,
     * And the second menu item terminates the application
     * */
    connect(viewWindow, SIGNAL(triggered()), this, SLOT(show()));
    connect(quitAction, SIGNAL(triggered()), this, SLOT(pad_exit()));

    menu->addAction(viewWindow);
    menu->addAction(quitAction);

    /* Set the context menu on the icon
     * And show the application icon in the system tray
     * */
    trayIcon->setContextMenu(menu);
    trayIcon->show();

    /* Also connect clicking on the icon to the signal processor of this press
     * */
    connect(trayIcon, SIGNAL(activated(QSystemTrayIcon::ActivationReason)),
            this, SLOT(iconActivated(QSystemTrayIcon::ActivationReason)));

}

void pad::iconActivated(QSystemTrayIcon::ActivationReason reason)
{
    switch (reason){
    case QSystemTrayIcon::Trigger:
        /* The event is ignored if the checkbox is not checked
         * */
        if(!this->isVisible()){
            this->show();
        } else {
            this->hide();
        }
        break;
    default:
        break;
    }
}

void pad::pad_show()
{

}

void pad::pad_kill()
{
    QMessageBox::warning(this, "Внимание","Программа не зарегистрирована");
    exit(1);
}

void pad::pad_key()
{
    this->setEnabled(false);
//    this->hide();
    m_key = new QDialog(0,0);
    m_ui_check=new Ui::check;
    m_ui_check->setupUi(m_key);
    QObject::connect(m_ui_check->pb_ok, SIGNAL(clicked()),this,SLOT(pad_key_enter()));
    QObject::connect(m_ui_check->pb_cancel, SIGNAL(clicked()),this,SLOT(pad_key_exit()));
    m_key->setWindowFlags(Qt::CustomizeWindowHint | Qt::WindowTitleHint);
    m_key->setFixedSize(m_key->geometry().width(),m_key->geometry().height());
    //    m_about->setEnabled(false);
    m_key->show();
}

void pad::pad_minimize()
{
    this->hide();
}

void pad::pad_tutorial()
{
    QString docfile=QCoreApplication::applicationDirPath() + "/tutorial.html";
    if (!QDesktopServices::openUrl("file:///" + docfile))
        QMessageBox::warning(this, "Ошибка", "Не могу открыть браузер");
}

void pad::closeEvent(QCloseEvent * event)
{
    pad_exit();
}
void pad::pad_exit()
{
    if(m_dirty)
    {
        QMessageBox msgBox;  //www.itmathrepetitor.ru
        msgBox.setText("Информация");
        msgBox.setInformativeText("Проект не сохранен. Сохранить?");
        msgBox.setStandardButtons(QMessageBox::Ok | QMessageBox::Cancel);
        msgBox.setIcon(QMessageBox::Information);
        msgBox.setDefaultButton(QMessageBox::Ok);
        int res = msgBox.exec();
        if (res == QMessageBox::Ok) //нажата кнопка Ok
            pad_save();
    }
    exit(0);
}

void pad::pad_solve()
{
    this->setEnabled(false);
//    this->hide();
    m_solve = new QDialog(0,0);
    m_ui_solve=new Ui::solve;
    m_ui_solve->setupUi(m_solve);
    QObject::connect(m_ui_solve->pb_find, SIGNAL(clicked()),this,SLOT(pad_solve_do()));
    QObject::connect(m_ui_solve->pb_cancel, SIGNAL(clicked()),this,SLOT(pad_solve_exit()));
    m_solve->setWindowFlags(Qt::CustomizeWindowHint | Qt::WindowTitleHint);
    m_solve->setFixedSize(m_solve->geometry().width(),m_solve->geometry().height());
    //    m_about->setEnabled(false);
    m_solve->show();

}

void pad::pad_solve_do(){
    QString fx=m_ui_solve->le_fx->text();
    double a,b,rv;
    bool ok;
    a=m_ui_solve->le_a->text().toDouble(&ok);
    if(!ok){
        return;
    }
    b=m_ui_solve->le_b->text().toDouble(&ok);
    if(!ok){
        return;
    }
    //    func_params p;
    m_ex=fx.toStdString();
    auto func=[ this ] ( double x,void * sd)
    {
        double ev;
        //        printf("ex=%s\n",m_ex.c_str());
        te_expr f;
        std::vector<te_variable> vars={{"x", &x}};
        te_parser parser;
        parser.set_vars(vars);
        parser.compile(m_ex.c_str());
        return parser.evaluate();
    };
    double sep=(b-a)/1000;
    QString r_text;
    std::vector<double> roots;
    double x_c;
    int i=0,nf=0;
    int cnt=1;
    x_c=FindNZero(func,a,b,NULL,sep,i,&nf);
    if(!nf){
        roots.push_back(x_c);
        r_text=r_text+QString("x[")+QString::number(cnt)+QString("]=")+QString::number(x_c,'g',14)+QString("\n");
        cnt++;
    }
    nf=0;
    i++;
    while(!nf){
        x_c=FindNZero(func,a,b,NULL,sep,i,&nf);
        if(!nf){
            roots.push_back(x_c);
            r_text=r_text+QString("x[")+QString::number(cnt)+QString("]=")+QString::number(x_c,'g',14)+QString("\n");
            cnt++;
        }
        i++;
    }
    m_ui_solve->te_roots->setText(r_text);
}

void pad::pad_solve_exit(){
    m_solve->close();
//    this->show();
    this->setEnabled(true);
}

void pad::pad_plot_exit()
{
    m_plot->close();
//    this->show();
    this->setEnabled(true);
}

void pad::pad_plot_export()
{
    auto customPlot=m_ui_plot->view;
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save plot"), QString(), tr("Graph Files (*.jpg)"));
    if(fileName.isEmpty()){
        return;
    }
    customPlot->saveJpg(fileName);
}

bool pad::pad_wait()
{
#if defined Q_OS_LINUX
#elif defined Q_OS_WIN
#endif
}

void pad::pad_load()
{
    char str[2048];
    pad_clear();
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Session"), QString(), tr("Session Files (*.dat *.log *.txt)"));
    if(fileName.isEmpty()){
        return;
    }
    FILE *f;
    f=fopen(fileName.toStdString().c_str(),"r");
    while(fgets(str,2000,f)!=NULL)
    {
        if(str[strlen(str)-1]=='\n'){
            str[strlen(str)-1]='\0';
        }
        m_expression=std::string(str);
        evaluate();
    }
    fclose(f);
}

void pad::pad_tray(){

}

void pad::pad_save()
{
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save session"), QString(), tr("Session Files (*.dat *.log *.txt)"));
    if(fileName.isEmpty()){
        return;
    }
    FILE *f;
    f=fopen(fileName.toStdString().c_str(),"w");
    for(int i=0;i<m_history.size();i++){
        fprintf(f,"%s\n",m_history[i].c_str());
    }
    fclose(f);
}

void pad::pad_copy()
{
    ui->plainTextEdit->copy();
}

void pad::pad_plot_do(){
    //    char *ex;
    QString fx=m_ui_plot->le_fx->text();
    double xmin,xmax,ymin,ymax,rv;
    bool ok;
    xmin=m_ui_plot->le_xmin->text().toDouble(&ok);
    if(!ok){
        return;
    }
    xmax=m_ui_plot->le_xmax->text().toDouble(&ok);
    if(!ok){
        return;
    }
    ymin=m_ui_plot->le_ymin->text().toDouble(&ok);
    if(!ok){
        return;
    }
    ymax=m_ui_plot->le_ymax->text().toDouble(&ok);
    if(!ok){
        return;
    }
    //    func_params p;
    m_ex=fx.toStdString();
    //    printf("a=%lg b=%lg expr=%s\n",a,b,m_ex.c_str());
    auto func=[ this ] ( double *x,void * sd)
    {
        double ev;
        //        printf("ex=%s\n",m_ex.c_str());
        te_expr f;
        std::vector<te_variable> vars={{"x", &x[0]}};
        te_parser parser;
        parser.set_vars(vars);
        parser.compile(m_ex.c_str());
        return parser.evaluate();
    };
#define SIZE_PT 1000
    QVector<double> x, y; // initialize with entries 0..100
    double x_t,y_t;
    for (int i=0; i<SIZE_PT; ++i)
    {
        x_t = xmin+i*(xmax-xmin)/(SIZE_PT-1);
        ok=true;
        try{
            y_t = func(&x_t,nullptr);
        }
        catch(...)
        {
            ok=false;
        }
        if(ok)
        {
            x.push_back(x_t);
            y.push_back(y_t);
        }
    }
    auto customPlot=m_ui_plot->view;
    // create graph and assign data to it:
    customPlot->addGraph();
    customPlot->graph(0)->setData(x, y);
    // give the axes some labels:
    customPlot->xAxis->setLabel("x");
    customPlot->yAxis->setLabel("y");
    // set axes ranges, so we see all data:
    customPlot->xAxis->setRange(xmin, xmax);
    customPlot->yAxis->setRange(ymin, ymax);
    customPlot->replot();
}

void pad::pad_plot()
{
//    this->hide();
    this->setEnabled(false);
    m_plot = new QDialog(0,0);
    m_ui_plot=new Ui::plot;
    m_ui_plot->setupUi(m_plot);
    QObject::connect(m_ui_plot->pb_plot, SIGNAL(clicked()),this,SLOT(pad_plot_do()));
    QObject::connect(m_ui_plot->pb_cancel, SIGNAL(clicked()),this,SLOT(pad_plot_exit()));
    QObject::connect(m_ui_plot->pb_export, SIGNAL(clicked()),this,SLOT(pad_plot_export()));
    m_plot->setWindowFlags(Qt::CustomizeWindowHint | Qt::WindowTitleHint);
    m_plot->setFixedSize(m_plot->geometry().width(),m_plot->geometry().height());
    //    m_about->setEnabled(false);
    m_plot->show();

}

void pad::pad_paste()
{
    ui->lineEdit->paste();
}

void pad::pad_copylast()
{
    setSelection(QString::number(m_ans,'g',14));
}

#ifdef Q_OS_WINDOWS
bool pad::nativeEvent(const QByteArray &eventType, void *message, long *result)
{
    Q_UNUSED(eventType)
    Q_UNUSED(result)
    MSG* msg = reinterpret_cast<MSG*>(message);
    if(msg->message == WM_HOTKEY){
        if(msg->wParam == 100){
            clip_scan();
            return true;
        }
    }
    return false;
}
#endif

void pad::pad_about()
{
//    this->hide();
    this->setEnabled(false);
    m_about = new QDialog(0,0);
    auto ui=new Ui::about;
    ui->setupUi(m_about);
    QObject::connect(ui->pushButton, SIGNAL(clicked()),this,SLOT(pad_about_exit()));
    m_about->setWindowFlags(Qt::CustomizeWindowHint | Qt::WindowTitleHint);
    m_about->setFixedSize(m_about->geometry().width(),m_about->geometry().height());
    //    m_about->setEnabled(false);
    m_about->show();
}

void pad::pad_integrate()
{
//    this->hide();
    this->setEnabled(false);
    m_integrate = new QDialog(0,0);
    m_ui_int=new Ui::integrate;
    m_ui_int->setupUi(m_integrate);
    QObject::connect(m_ui_int->pb_eval, SIGNAL(clicked()),this,SLOT(pad_eval_integrate()));
    QObject::connect(m_ui_int->pb_cancel, SIGNAL(clicked()),this,SLOT(pad_integrate_exit()));
    //    m_about->setWindowFlags(Qt::CustomizeWindowHint | Qt::WindowTitleHint);
    //    m_about->setFixedSize(m_about->geometry().width(),m_about->geometry().height());
    //    m_about->setEnabled(false);
    m_integrate->show();
}

QString pad::selectedText()
{
    QString selectedText;
#if defined(Q_OS_LINUX)
    selectedText = QApplication::clipboard()->text(QClipboard::Selection);
#elif defined(Q_OS_WIN) // Send Ctrl + C to get selected text
    // Save original clipboard data
    QThread::msleep(500);
    QVariant originalClipboard;
    if (QApplication::clipboard()->mimeData()->hasImage())
        originalClipboard = QApplication::clipboard()->image();
    else
        originalClipboard = QApplication::clipboard()->text();

    // Generate key sequence
    INPUT copyText[4];

    // Set the press of the "Ctrl" key
    copyText[0].ki.wVk = VK_CONTROL;
    copyText[0].ki.dwFlags = 0; // 0 for key press
    copyText[0].type = INPUT_KEYBOARD;

    // Set the press of the "C" key
    copyText[1].ki.wVk = 'C';
    copyText[1].ki.dwFlags = 0;
    copyText[1].type = INPUT_KEYBOARD;

    // Set the release of the "C" key
    copyText[2].ki.wVk = 'C';
    copyText[2].ki.dwFlags = KEYEVENTF_KEYUP;
    copyText[2].type = INPUT_KEYBOARD;

    // Set the release of the "Ctrl" key
    copyText[3].ki.wVk = VK_CONTROL;
    copyText[3].ki.dwFlags = KEYEVENTF_KEYUP;
    copyText[3].type = INPUT_KEYBOARD;

    // Send key sequence to system
    SendInput(4, copyText, sizeof(INPUT));

    // Wait for the clipboard to change
    QEventLoop loop;
    QTimer timer; // Add a timer for the case where the text is not selected
    loop.connect(QApplication::clipboard(), &QClipboard::changed, &loop, &QEventLoop::quit);
    loop.connect(&timer, &QTimer::timeout, &loop, &QEventLoop::quit);
    timer.start(1000);
    loop.exec();

    // Translate the text from the clipboard if the selected text was not copied
    if (timer.isActive())
        return QApplication::clipboard()->text();
    else
        timer.stop();

    // Get clipboard data
    selectedText = QApplication::clipboard()->text();

    // Return original clipboard
    if (originalClipboard.type() == QVariant::Image)
        QApplication::clipboard()->setImage(originalClipboard.value<QImage>());
    else
        QApplication::clipboard()->setText(originalClipboard.toString());
#endif
    return selectedText;
}

void pad::setSelection(QString text)
{
    QClipboard* clipboard = QApplication::clipboard();

    clipboard->setText(text, QClipboard::Clipboard);

    if (clipboard->supportsSelection()) {
        clipboard->setText(text, QClipboard::Selection);
    }
    QThread::msleep(750); //workaround for copied text not being available...
}



void pad::clip_scan()
{
    const char *str;
    int len;
    QString text=selectedText();
    QString result;
    m_expression=text.toStdString();
    evaluate();
    result=QString::number(m_ans,'g',14);
    setSelection(result);
#if defined Q_OS_LINUX
    nativeEventFilter->emu_out();
#elif defined Q_OS_WIN
    INPUT copyText[4];

    // Set the press of the "Ctrl" key
    copyText[0].ki.wVk = VK_CONTROL;
    copyText[0].ki.dwFlags = 0; // 0 for key press
    copyText[0].type = INPUT_KEYBOARD;

    // Set the press of the "C" key
    copyText[1].ki.wVk = 'V';
    copyText[1].ki.dwFlags = 0;
    copyText[1].type = INPUT_KEYBOARD;

    // Set the release of the "C" key
    copyText[2].ki.wVk = 'V';
    copyText[2].ki.dwFlags = KEYEVENTF_KEYUP;
    copyText[2].type = INPUT_KEYBOARD;

    // Set the release of the "Ctrl" key
    copyText[3].ki.wVk = VK_CONTROL;
    copyText[3].ki.dwFlags = KEYEVENTF_KEYUP;
    copyText[3].type = INPUT_KEYBOARD;

    // Send key sequence to system
    SendInput(4, copyText, sizeof(INPUT));
#endif
}

void pad::pad_about_exit()
{
    m_about->close();
//    this->show();
    this->setEnabled(true);
}
void pad::pad_key_exit()
{
    m_key->close();
//    this->show();
    this->setEnabled(true);
}

void pad::pad_check()
{
    QSettings registry(QSettings::NativeFormat, // Format
                       QSettings::UserScope, // Scope
                       "IOFFE", // Organization
                       "SmallPAD" // Application
                       );
    QStringList keyList = registry.allKeys();
    QStringList valList;
    for(int i=0;i<keyList.size(); i++)
        valList.push_back(registry.value(keyList[i]).toString());
    if(keyList.size()==0)
    {
        if(!m_t->isActive())
        {
            m_t->start(300000);
        }
    }
    else
    {
        char tmp[512];        
        QString name = keyList.at(0);
        QString pass = valList.at(0);
        uint32_t crc = crc32((unsigned char *)name.toStdString().c_str(),(size_t)strlen(name.toStdString().c_str()));
        double x=(double) crc/(double)0xFFFFFFFF;
        for(int i=0;i<100;i++)
        {
            x=4*x*(1-x);
        }
        sprintf(tmp,"%.9lf",x);
        if((strlen((char *)pass.toStdString().c_str())!=9) || (strlen(name.toStdString().c_str())==0)){
            if(!m_t->isActive())
            {
                m_t->start(300000);
            }
        }

        bool test = true;
        for(int i=0;i<9;i++){
            if(((char *)pass.toStdString().c_str())[i]!=tmp[i+2])
            {
                test = false;
                if(!m_t->isActive())
                {
                    m_t->start(300000);
                }
            }
        }
        if(test)
        {
            m_t->stop();
        }
    }
}

void pad::pad_key_enter()
{
    QSettings registry(QSettings::NativeFormat, // Format
                       QSettings::UserScope, // Scope
                       "IOFFE", // Organization
                       "SmallPAD" // Application
                       );
    registry.clear();
    QString name = m_ui_check->le_user->text();
    QString pass = m_ui_check->le_pass->text();
    registry.setValue(name,pass);
    pad_check();
    m_key->close();
//    this->show();
    this->setEnabled(true);
}

void pad::pad_integrate_exit()
{
    m_integrate->close();
    this->setEnabled(true);
//    this->show();
}

void pad::pad_eval_integrate()
{
    QString fx=m_ui_int->le_fx->text();
    double a,b,rv;
    bool ok;
    a=m_ui_int->le_left->text().toDouble(&ok);
    if(!ok){
        return;
    }
    b=m_ui_int->le_right->text().toDouble(&ok);
    if(!ok){
        return;
    }
    //    func_params p;
    m_ex=fx.toStdString();
    //    printf("a=%lg b=%lg expr=%s\n",a,b,m_ex.c_str());
    auto func=[ this ] ( double *x,void * sd)
    {
        double ev;
        //        printf("ex=%s\n",m_ex.c_str());
        te_expr f;
        std::vector<te_variable> vars={{"x", &x[0]}};
        te_parser parser;
        parser.set_vars(vars);
        parser.compile(m_ex.c_str());
        return parser.evaluate();
    };
    rv=GaussIntegrate(func,nullptr,1,&a,&b,10000);
    m_ui_int->le_result->setText(QString::number(rv,'g',14));
}

bool pad::eventFilter(QObject* obj, QEvent *event)
{
    if (obj == ui->lineEdit)
    {
        if (event->type() == QEvent::KeyPress)
        {
            QKeyEvent* keyEvent = static_cast<QKeyEvent*>(event);
            if (keyEvent->key() == Qt::Key_Up)
            {
                if(m_ptr==-1){
                    m_ptr=m_history.size()-1;
                }else if(m_ptr==0){
                    m_ptr=0;
                }else{
                    m_ptr--;
                }
                if(m_ptr==-1){
                    ui->lineEdit->setText(QString(""));
                }else{
                    ui->lineEdit->setText(QString::fromStdString(m_history[m_ptr]));
                }
                return true;
            }
            else if(keyEvent->key() == Qt::Key_Down)
            {
                if(m_ptr==-1){
                    m_ptr=-1;
                }else if(m_ptr==m_history.size()-1){
                    m_ptr=-1;
                }else{
                    m_ptr++;
                }
                if(m_ptr==-1){
                    ui->lineEdit->setText(QString(""));
                }else{
                    ui->lineEdit->setText(QString::fromStdString(m_history[m_ptr]));
                }
                return true;
            }
        }
        return false;
    }
    return QWidget::eventFilter(obj, event);
}

void pad::pad_clear()
{
    ui->plainTextEdit->clear();
    m_history.clear();
    count=0;
    m_dirty = false;
}

pad::~pad()
{
    delete ui;
}

void pad::evaluate()
{
    m_dirty = true;
    //    te_parser parser;
    te_expr expr;
    int err,i,j;
    double result;
    std::string var;
    QString output;
    m_vars.clear();
    m_vars.push_back({"ans", &m_ans});
    //Check for variables;
    char *left,*right;
    left=(char *)malloc(m_expression.size()*sizeof(char));
    right=(char *)malloc(m_expression.size()*sizeof(char));
    bool flag=false;
    for(i=0;i<m_expression.size();i++){
        if(m_expression[i]=='='){
            flag=true;
        }
    }
    //Do split;
    if(flag){
        i=0;
        while(m_expression[i]!='=')
        {
            left[i]=m_expression[i];
            i++;
        }
        left[i]='\0';
        i++;
        for(j=0;i<m_expression.size();i++,j++){
            right[j]=m_expression[i];
        }
        right[j]='\0';
    }
    bool exist=false;
    int idx;
    if(flag){
        var=std::string(left);
        //Check for existence;
        for(i=0;i<m_vars_string.size();i++){
            if(var.compare(m_vars_string[i])==0){
                exist=true;
                idx=i;
            }
        }
        if(!exist){
            m_vars_string.push_back(var);
            m_vars_double.push_back(0.0);
        }
    }
    for(i=0;i<m_vars_string.size();i++){
        m_vars.push_back({m_vars_string[i].c_str(),&m_vars_double.data()[i]});
    }
    parser.set_vars(m_vars);
    if(!flag){
        parser.compile(m_expression.c_str());
    }else{
        parser.compile(right);
    }
    result=parser.evaluate();
    m_ans=result;
    if(flag){
        if(!exist){
            m_vars_double[m_vars_double.size()-1]=result;
        }else{
            m_vars_double[idx]=result;
        }
    }
    m_results.push_back(result);
    m_history.push_back(m_expression);
    count++;
    output=QString("In[")+QString::number(count)+QString("]=")+QString::fromStdString(m_expression);
    ui->plainTextEdit->appendPlainText(output);
    output=QString("Out[")+QString::number(count)+QString("]=")+QString::number(result,'g',15);
    ui->plainTextEdit->appendPlainText(output);
}


///**
// * @brief pad::eval_func
// * Производит вычисление функции по строке.
// * @param x
// * @param sd
// * @return
// */
//static double pad::eval_func(double *x,void *sd)
//{
//    func_params *ptr;
//    ptr=(func_params *)sd;
//    te_expr f;
//    std::vector<te_variable> vars;
//    m_vars.push_back({"x", &x[0]});
//    te_parser parser;
//    parser.set_vars(vars);
//    parser.compile(ptr->expression);
//    return parser.evaluate();
//}

void pad::on_lineEdit_returnPressed()
{
    m_expression=ui->lineEdit->text().toStdString();
    evaluate();
    m_ptr=-1;
    ui->lineEdit->clear();
}
