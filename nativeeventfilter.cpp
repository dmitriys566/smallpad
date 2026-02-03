#include "nativeeventfilter.h"
#include <QVector>
#include <QX11Info>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <xcb/xcb.h>
#include <X11/keysym.h>
#include <X11/extensions/XTest.h>
#include <X11/Xatom.h>
//#include <X11/Xmu/Atoms.h>
//#include <chrono>
//#include <thread>

extern "C" {
#include <xdo.h>
}

namespace
{
Display * m_display;        // Соединение с сервером X11
Window m_win;               // Захватываемое окно - в данном случае будет вся система
int keycode;                // код клавиши
unsigned int modifier;      // код модификаторов
/* Вектор дополнительных модификаторов Num Lock, Caps lock
     * Они тоже учитываются в X11, поэтому понадобяться все возможные комбинации
     * */
QVector<quint32> maskModifiers(){
    return QVector<quint32>() << 0 << Mod2Mask << LockMask << (Mod2Mask | LockMask);
}
}

NativeEventFilter::NativeEventFilter(QObject *parent) : QObject(parent)
{
    m_display = QX11Info::display();        // Создадим подключение к серверу
    m_win = DefaultRootWindow(m_display);   // и вытащим из него захватываемое окно с помощью макроса
}

bool NativeEventFilter::nativeEventFilter(const QByteArray &eventType, void *message, long *result)
{
    Q_UNUSED(eventType)
    Q_UNUSED(result)

    /* В вот обработка события строится уже на библиотеке XCB вместо Xlib.
         * Вроде как, получая событие Qt знает его в качестве XCB события,
         * но не знает его в качестве события Xlib, хотя использовать более
         * легкий синтаксис Xlib для установки хоткеев нам никто не запрещает
         * */
    xcb_key_press_event_t *keyEvent = 0;

    // Итак проверяем, что это было xcb событие
    if (eventType == "xcb_generic_event_t") {
        // кастуем сообщение собственно в xcb событие
        xcb_generic_event_t *event = static_cast<xcb_generic_event_t *>(message);

        // проверяем, что произошло нажатие клавиши
        if ((event->response_type & 127) == XCB_KEY_PRESS){

            // Если так, то кастуем сообщение в событие нажатия клавиши
            keyEvent = static_cast<xcb_key_press_event_t *>(message);

            // Далее проверям, является ли это событие нужным хоткее
            foreach (quint32 maskMods, maskModifiers()) {
                if((keyEvent->state == (modifier | maskMods ))
                        &&  keyEvent->detail == keycode){
                    emit activated();   // и посылаем сигнал
                    return true;
                }
            }
        }
    }
    return false;
}

void NativeEventFilter::setShortcut()
{
    unsetShortcut();        /* Вначале для перестраховки отключим предполагаемый хоткей,
                                 * даже несмотря на то, что будет мусор в первый раз в параметрах хоткея.
                                 * */

    // получим код клавиши из KeySym определения и соединения с сервером X11
    keycode = XKeysymToKeycode(m_display, XK_equal);
    modifier = ControlMask | Mod1Mask; // Зададим модификатор

    /* А теперь пройдемся по всем возможным комбинациям с учётом Num Lock и Caps Lock
         * устанавливая хоткеи
         * */
    foreach (quint32 maskMods, maskModifiers()) {
        XGrabKey(m_display,         // указываем соединение с X11
                 keycode ,          // код клавиши
                 modifier | maskMods,   // модификатор со всеми возможными масками
                 m_win,             // Захватываемое окно
                 True,              // Является ли приложение владельцем события. в данном примере не принципиально.
                 GrabModeAsync,     // Обязательно Ассинхронный режим обработки, иначе, рискуете встрять
                 GrabModeAsync);    // с замороженной системой, не реагирующей ни на какие воздействия, если
        // заранее не напишите корректную передачу события обратно в систему,
        // а скорее всего так и будет
    }
}

void NativeEventFilter::unsetShortcut()
{
    // Проходим по всем возможным комбинациям и удаляем хоткей
    foreach (quint32 maskMods, maskModifiers()) {
        XUngrabKey(m_display,
                   keycode,
                   modifier | maskMods,
                   m_win);
    }
}

//void NativeEventFilter::emu_in()
//{
//    xdo_t * x = xdo_new(NULL);
//    sleep(1);
//    xdo_send_keysequence_window(x, CURRENTWINDOW, "Ctrl+x", 200000);
//    xdo_free(x);
//}

void NativeEventFilter::emu_out()
{
        xdo_t * x = xdo_new(NULL);
        sleep(1);
        xdo_send_keysequence_window(x, CURRENTWINDOW, "Ctrl+v", 200000);
        xdo_free(x);
}

bool NativeEventFilter::wait(){
    char keys_return[32];
    bool value=false;
    XQueryKeymap( m_display, keys_return );
    for(int i=0;i<32;i++){
        if (keys_return[i]!=0){
            value=true;
        }
    }
    return value;
}

//void NativeEventFilter::setWindowText(QString text){

//}

