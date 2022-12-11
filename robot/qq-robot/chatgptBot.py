import flask, json
from flask import request
from revChatGPT.revChatGPT import Chatbot

config = {
    "session_token":
    "eyJhbGciOiJkaXIiLCJlbmMiOiJBMjU2R0NNIn0..ek_xGOQV6_S8pky_.vYC-b8UbX5u6mOzF6F_e-eCTD7eJkprHEK1BEQ0QGskVymAkvhRh-Eogt-APdneuu2zxTPegxgKUSYYAri5SxWIciVoldeO6n2WNkrofM0OzgLKTHfc-w3aateBBthp0aspvB5-Spz8xgPuiI_-V5AAUvI92cjPxcH5mb9B-RPVXeNhBCAsaS7NiMe4zrUwK_9wiXbf5DPwHHyTg1msXaOT_sfCfUyUfbrk1lsffa0gJ4wcLev8N__WtwdtbZKnwALHoORfrCBxo83QAqWgP22fNIWASaNZf5brOvE04fZwdrXI55pJKssVfmxhtFA4aJYZRNYnFc8wU7d0mV9_UqD6E4N6R3nLmx1efd6Nl-C796WVdxkRfnq7ZbjGuetiSdAf-n9L0XNJOHSWp1CyTAeGmXSVI7V7c41cbs3CqRUEjvCmNA9hSseUZrtq4qKXSm-wZMgf_UxUS33Z2ygmlxUjj_1-xFn8dSIJ0ZpQ1lxQfeqKUb1Iar4Hwk4JDa_cihjYsKXrRQdjxGpPcBVVEDbZNXGx8BshZDGmlzwUdLhL7TDfZf8C9BWpx-Daj-ZXmRPEdZouwMrsCDF5c_StPYFkDhnUd7HtWICG1-zb3dS1rT2Lvt5BkZXp7Cs2jIUUfHTrDVXH20PDJt2aGQ4HI5i-P5xj4NObwmPJUxucgIxTtve0FaYKrbeg9euBTIcolmXORhgtXGRdQTg_h5Pa2B7Kik1FGb2a58UctdRHkvanNVL96-SiGMY7LswNehO1qVh5skzCZaS2YtOMRhhJoQxenPLyOzf7H4DD-mEBUiuwqR15VCXkgs2BCbNuxlBIcT9S6D4lektPVsONvReikCY2OMt8K15DYMBI7n8rybUO2cSRLng7zY2j30ZI35QWownVC4KvpBBF2dShr9JSKIu5f_XMDLx6ak3SrSTeBN8qUUN8vb8X42Sa4cz2a1SGfgCBh6Vs6dTB_K1H87_QcM05gLdoMk98Kj_P1IOUR3xEtIo2L-Ff2raUW5Szpd4MxOI8fZiTF2gecHXQjDBn82sAKSWRalTSR7QZHLi2nO0J4svpuuQqnKsKQddcc3lr50PJUyc1a9OsHDY0HhhhJGmwD-GcC-PC-AMQtwGMjfjsFrXrjmu3ekUwy2xMLNFiA0pYOg_QU03qckDPAL9Ro4nSk--Hk_v0HZC4rFgxlzmM5TtjRR4So53vkL_8VanyAAq1LkPtmIklM8-jRCHjuezEgQtm-TO-QkoNRmDeoXIQwEJ9JsQxFDV5ZcKpK4w4m9Zocs-hyn04JOKfkpzzTDkkzEhegKkaZ3LpZZJP_KJvZoWoezekw5JXipqa3CDmlGoGrfcLRaC9z-WTxghCnlSvkMNIRfI8y3PjqgSHBLo6P_sFWZHUqTEhb-h4JGzbin8ATGEUt5P1QS5t-s_lKoT7mtXBu78-lFDaZ0BxGoGTA-EE6h2PcVtxm0R-_KPw2dIUG5w-3UoGpRan3ARFeI52LboNqKZVC61g62n7hdgF2SmJ--qr71dSkJU0yzjdUQS6GzuYsuF71-6Lk2ApteCxE_bm_s6N5KQWqC-mrka-H8KxoG4lWx_jZ0udCBQZ0A6HQPTQmBGxkXI-a18dQWrkWJAYfm0eUDUG1fi3H0NZEEmYKWCdTnGBdoIC6DngDDhVjWRfCiUu2DRTKwYEYDO1OJb2KnYt-5MiwxdaFZoAFJPDSlm6b56wg44X_yR0GQfPnwNYa3QX6ujhIm6zkrbUctq06YYSguVMqDHKB7HxPeAJ85wo8yeUgpWNKwOTpFoxvtggcx5NJHZe7QF1vl6rMZOoYr9TVOPAjDap1-yOsjM_FkQbj62ssfjmB6sH689DKJqSWdOzlWXUt27MvB8IVFQqHrYG8f6genCmuAv_pi-VzGNCAtbUkW1DuwhaDnD5GfgB95M2r0AnWcNH2SPWmOHnKtv1Q5pc07HglaxywAvga0q1C0vi3RJ5gVyV_voOVlgdknobUr1Xdc5CXR_xp1VDjhM7fwjgQYXRHrcTpKCvlOQqLeVxagA8pAFIUSFchoty3VSnWW4jiDs3vCzTiubcPJGsoc1ILBRnmI8uROrjea-KnvrZl0JE1Cuy2ECFQItgHU-Y0B-Cw_k8dQBHcyehQdspYLqMgGsCB6n6SrzDNqfcbrTn4UqdzyriN1UMnDC0ioOrsLNyzXk03aznybB6Lf8ebfhq0HTE6whTfEDSQ1G94piqC4yHxBbnrRRsFDv1qxn8a7kd5k8fGmclg5gdry6FiCdV8x52MB3CT9lYMkupLKSsYGlfw2vVPgFj1vD1XkO57AabtHL6DOwlV5YpbZA.u6jXMk-KSDYonbAJSgs8yQ"
}
# 创建一个服务，把当前这个python文件当做一个服务
server = flask.Flask(__name__)
chatbot = Chatbot(config, conversation_id=None)


def chat(msg):
    message = chatbot.get_chat_response(msg)['message']
    print(message)
    return message


@server.route('/chat', methods=['post'])
def chatapi():
    requestJson = request.get_data()
    if requestJson is None or requestJson == "" or requestJson == {}:
        resu = {'code': 1, 'msg': '请求内容不能为空'}
        return json.dumps(resu, ensure_ascii=False)
    data = json.loads(requestJson)
    print(data)
    try:
        msg = chat(data['msg'])
    except Exception as error:
        print("接口报错")
        resu = {'code': 1, 'msg': '请求异常: ' + str(error)}
        return json.dumps(resu, ensure_ascii=False)
    else:
        resu = {'code': 0, 'data': msg}
        return json.dumps(resu, ensure_ascii=False)


if __name__ == '__main__':
    server.run(port=7777, host='0.0.0.0')