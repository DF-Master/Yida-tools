import flask, json
from flask import request
from revChatGPT.revChatGPT import Chatbot

config = {
    "session_token":
    "eyJhbGciOiJkaXIiLCJlbmMiOiJBMjU2R0NNIn0..EQNgt9AEqF-Qd4nZ.Hx_LQM3MWe37_A_R4YBVJxmdQmrPLYhmt8amzllmXEB87hsORBt6HyPS8wlJlhELlO1V0XfGkjLWeUZD6bfv0bH93L0H1LgnZdMiyA7VBk8ZD_AzMrs9C_YXzU-fbXZQG6bNxA4L8CJau_euv0UIUNI07BVhdAHs254vQea0ri3Fm9U62fgAvRei9k1sRo4oDXoc-srk7DZ8mB64vn0QslSZET0prEnKsBB4tA081WzGkXC0Swex8XHhTDpnQNJeb-tdQb9A6tNV8gin43UBkm7f7hCH-aCG4rEo4d31eM1gARUHD7cvGzS0WpPMTg4kLpzIO5NPWHmllDN0aoDpGDvyRVNqumUR4dGcqB7cRTfvwMF1U2d89mNQ3rRoTUEx5qJAo76GKOQD2i-0COxuRKvdKKjpOyoJvr3Za8uO4ivekYPi_OKE2R-vQnf_NPP8nGS2WBV2Je3tdBfnFDx5WNwlbrTKWHl5iXjPsn8dNOJWTb_dqDvicJiSfPYQMk0AXWDfUovMBdAXwjXNaMhCGbViO4rXlWzuEz4OZtk3i399tpaH0ZZk1Y3MuiYU_peNi-_RREeZ1L-qTzyGIjjeUhlA42pyNNx2WG_fi5QFXGlyPbnDMxSGq8e2AP4Ok-eVn3qT8wkJ9EsmPnQhI9QkKufCEWAqua2iZAYxH4071aXK0_0N76ILJsoC23ho3HAhZ8kgwvrZExJo55lycqMS5evwdEx1NaDHOmiYTvE7HFfHFioTksnmI92xZJB--pptHuPribf756EhzRwzoAkfjg5xzdfupWG8xwBfHSk2ICMCRli8YitxkRVWpP7zcjkwR6_DqUhTfDxJxiz2YVwRATyVMSsq_D9-ejuDSlRi6UvSVNDNd9AmJG9hhlngI9szdYXNRhLsM9bvBYA9U-BJ1nkuv7Wn0ybVCsnu3xLFG1aE_ahVN_O3aCVnU3cJyk-JoZO0tSPTzGaBIkbmVB1hvyVlsIJDD0988XDvMjbyOa4QUUPTla1Eslg9mvL2EXPaVEKda4rpt3iZlGBlB7vaTf01AJXDVQoh1_g_Uiz49obNVP1CSI9I7tkCsiBNhBuAnhpwQNWCXkLxB3WBcAFSHm4ZXookrUHDgpVBRw9qygSiyZbRmz9AmZjz0hOuWw_2EvmuKH5zc9d3N7BVZ0SUBpQeoCYXKJV9hq42TKmphLcYRzVt_VMWzKCsW6TQUl1p9sUCcfwKfiUzqHT2SfSWMcaGhu5itEzJJS1F94vG06WmWCMvZH_uoJE_8CFHD2R90Y503gTofdu-gePJZykq-vKKYeEdrEnNsYZPyYORHng9b7Q49UqnvNM-yWrNejlp4sg-sIjUCHrRJxj0w5D1tM9fugzrxe8ZyMkHu7HugCrh7rrmnb3IxlMHO95fqG2Z3iTD6_fN3pZcWnVHOlFsJJKTzdmVvvX4NW9zEUU7hH5lc5Qm85sgWC5UzqEryj6yEYRzlxzrMYEepCo0WrEGotESNPHyJras2unbTgIOx7H1XXqPaYlzCeWc_SowlTcxJ-zPHsQquaWORKunMetYs9HkVHr6ybpWajrFytPPpBrMZQCHcrN6E2H9n6D68gMassMLm-EmPuyUnuTVBlTA-xopNAim28IhiOFh7eI0NFiBeKoT5KyUu4M7URW6UteaJgyog7PFC98d8V5h5k04KTOqnpA-OugFuDZLfgZsF8kd45LgkLEAqG77XAE4HV3mJU1B7Q5EKT8teYB8nyW1UeGbfpOad7GfwRrX5UdW2IE1kWhy9XrkfzXKU1s7hi4K_ir2FcgyZHJ-d-qSsBVviJ2cEFV4-2w1Tixm4v9EH1Lw7dp-Mal1nq8tql4hAwHFpVQOD2HMJUAJCkdm3FwzZPn1_hz604LCFbkRFQeRa31OgollTm9WD_FXyD5W23lXGjf88kNpEtTbfhDlhgClMRANIkwTvRVy2sUckb_XW-SmaBI6K5u6AI0-W7K7ea3oTjNAxKElmzh_q_2vnsbfQ7I1dnIdhlNhb4aknT7plkq0Mnq4Oec9aD37oyLz1ZfgNmDXi4bsc8xz_A6AidZXJ2YZRH982ZUPqysv__ceEhYPQe2SMOc6zxfIpklNPQqIp5QHj2xFi-YbXEiFgxNTTRWS_bOPZWxuGcyUW-Mj12tMAPPyrJd20fHM6OK7PnMBp4QDhVKztYjS0b3Pyhwl4FCsmr9g92ryrraQ_cekMSxkDYIu80Tuz8b-tp137DQTnWGQq9mgkwUyC89u76ThPTBYK2o4C1tbXhEjif5wWILmL491c3pqwNtblfRD3ia9R_fr-kiS-dOT6tulsM6gzIioMdRuhg.ujU_90WwrqYB7HRBTxB8vQ"
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