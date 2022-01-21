from wxauto import *
import schedule
import time
import os

os.environ["http_proxy"] = "http://127.0.0.1:7890"
os.environ["https_proxy"] = "http://127.0.0.1:7890"

# Init
# 获取当前微信客户端
wx = WeChat()

# 获取会话列表
wx.GetSessionList()

# Def
name_list = [
    'mh1to9', '我是好人', 'Swiring', '烨子', 'Dracoʕ•ٹ•ʔ', 'ND是可燃的', '刷牙学长', '黄雨菲',
    '宝盖头一个诞辰的辰', '失眠', 'if', '雨中山', '蕴藉', '我自欲为江海客', '子不语', 'Clemmensen',
    '罗乌有', '。', '马亦然 Ferdinand', '洛融', '不羡仙', '二茂铁二胺', '潘高翔', 'Stephen Allen',
    'chen', 'Victorrrrr'
]

msg_auto = "------" + """大家好，我是本群的“提醒报备小助手”。您今天报备了吗？您汇报异动情况了吗?小助手提醒您：实验千万条，报备第一条；出门不报备，禁闭两行泪。""" + "------"

msg_functions = "目前实现的功能：自动化定时发放通知；自动化提醒未报备人员（需要完善）；功能查询；"

self_intro = "0w0蔻你吉瓦~这是一个自用报备提醒Bot,每天提醒四次报备及出入京情况确认，随开发者的心情迭代，并且可能会随着某一次的更新去世qaq"

msg_end = "##########"


# 建立统计函数
def Count_call():
    global name_list
    name_unsign = name_list
    mem_stat = {}
    wx.LoadMoreMessage()

    msgs = wx.GetAllMessage
    msgs = msgs[::-1]
    for msg in msgs:
        print(msg)
        if msg[0] == "DFMaster" and msg[1] == "$$$Start Here$$$":
            break
        elif msg[0] == """SYS""" or msg[0] == """DFMaster""" or msg[
                0] == """江意达""":
            continue
        else:

            try:
                name_unsign.remove(msg[0])
            except:
                print("Name not find:" + msg[0])
            else:
                try:
                    mem_stat[msg[1]] = mem_stat[msg[1]] + "," + msg[0]
                    print(mem_stat)
                except:
                    mem_stat[msg[1]] = msg[0]
                    print(mem_stat)
    msg_sign_situation = ""
    for type in mem_stat:
        msg_sign_situation += str(type) + ": " + str(mem_stat[type]) + ";;;"
        print(str(type) + ": " + str(mem_stat[type]))
    print(msg_sign_situation)
    print(name_list)
    wx.SendMsg(msg_sign_situation)
    wx.SendMsg("未报备人员:" + " ".join(name_unsign))


def Sign_Helper():

    # 向某人发送消息
    who = '【通知】21研一对一联系第二小组'

    wx.ChatWith(who)
    wx.SendMsg(msg_auto)
    Count_call()


def Start_Count():
    msg_block = "$$$Start Here$$$"
    wx.SendMsg(msg_auto)
    wx.SendMsg(msg_block)


def Listen_Order():
    msgs = wx.GetAllMessage

    for msg in msgs[::-1]:
        if msg[1] == msg_end:
            break
        elif msg[1] == "ID":
            wx.SendMsg(" ".join(name_list))
            wx.SendMsg(msg_end)
        elif msg[1] == "功能":
            wx.SendMsg(msg_functions)
            wx.SendMsg(msg_end)
        elif msg[1] == "介绍":
            wx.SendMsg(self_intro)
            wx.SendMsg(msg_end)
        elif msg[1] == "提醒":
            Count_call()
            wx.SendMsg(msg_end)
        else:
            continue


# 输出当前聊天窗口聊天消息
# msgs = wx.GetAllMessage
# print(msgs)
# for msg in msgs:
#     print('%s : %s' % (msg[0], msg[1]))

# Count_call()

if __name__ == '__main__':
    # 获取当前微信客户端
    wx = WeChat()

    # 获取会话列表
    wx.GetSessionList()

    # Count_call()
    schedule.every().day.at("00:00").do(Start_Count)

    schedule.every().day.at("09:30").do(Sign_Helper)
    schedule.every().day.at("12:00").do(Sign_Helper)
    schedule.every().day.at("16:50").do(Sign_Helper)
    # schedule.every(1).minutes.do(Listen_Order)

    while True:
        schedule.run_pending()  # 运行所有可以运行的任务
        time.sleep(30)