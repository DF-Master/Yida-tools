import flask, json
from flask import request
from revChatGPT.revChatGPT import Chatbot

config = {
    "session_token":
    "eyJhbGciOiJkaXIiLCJlbmMiOiJBMjU2R0NNIn0..Px-cthvqlUKI2UXq.CZo8gqFpetcaMOBbX6rvLdvpp2Kpajwqwt5bWFYLlgC_UuWxLEFFllHnJtI-3CDfpKLISQCpxHHGScDVw2_QQGwg_7NV9fN770z7bOd0zHo9Ymn_vdS9bOdrz6HpeHdl-2cgoKnnY10b3eoiBxWCdlGgZ-gQgq4Uvw69DA8_hi00hK-i9Sf-GvrY8BpfSdNdx92BMkO12S_KGVx_4u39bhcwQa1Es21kQQxJlmTbzpj4fkgRgW6rDWjD7fgF32_RNHXEZ0iZame8PtfEy1ylAkpuLBQuE6z4sCbE-5pdrx8RQyHwSTENSEUD8Z-hNJt3GKTW4Gyg-Jb8bhswuAUsPfG8kr0qJkB4mMBWnhBmCk-go42wqY4BAWQBRDFc0quZMt2ooyABU2RBaJgF5jtkUn-T2XVs6UZKBNmp6D3LBAbaXv0ycUjho16cqHo4O7vI_ijwdUYHwQURNdDnGISTeNGifMU0Q9nBD1-I8bfwRh5RDPL-WtOHsRwo4JKanTP2TbP8btKBEdcc2LNFHxLm4SV3yN959Fb8X5_uNJU5Or4nH2KOv01hSN74fxJclTnGd6165kgBgI3muEdQWVr65hJM7ddZ76QIZ4oHs3ZWPzyZq4HVb_aubD2yMM2oe_XXd3UjZ7penKF0v4EvaUi4NAhIxQddWbNppMBhfMQTZajAHj5moC9Ft-QRpryJqKYb6sf7CccMT-bjBQpyg9cBmqAdW7_Rcw79PCybuKdclhmcMbKOI0o3OPgWSEyEiQQCEpClVS9jlIsiGwcg6RdDAWR8um9_DwbsLV95aPNrC6PrcGUv1xUaO3NgvB56mEs55L8FOjh1ceWs6By9fKVIlaE62fWhUP_FKxTkJWwLWRn51_U9ItmoN2jxk1ieRMmV1GSjGQFzeXHeK08YxE1Yie4sMF_ocUlrvOzy0DJmQuRi4QOzDZKjs7TpFh9-yfUFb6yJANvrwejM_2jBjHoArqTvHvQEU9Al6ejCpkHy-vgxAG4PVzNp2D3Fq85BNGKJj0Aie_kaOCoBbfwzZ3MmaBnOCJuWOvaUB99otnZJrGHAu1YEEhBwYH5toTAOmvaRmgcQaKGzljMwXXjNBS8iO72e9CVmzVfPuMLhE8_5nsnX_SnkRV1Idfc6JZ7eEgFCHubBlMn7KDJtx25Pa72AECLQZpegPrq0WErVxor2OYofsuTcmnk_ZEYvhUfKBRu0trB0aOiMoRZAvoaniMv86qq-ra9RORhcaNrFLZOOHB5CdIfrX9Gi0BVJls697t-xE2Mumz8bu5DrIE7aAIP5HSQzORF611P2tcs9oeQaxvbQyOZyV2Dt5NtBvGGSFn4x0pRaZvnQ90MB14EHcfTj9ynUw9CW5wb1xVATTHWbpcQL-tSTtRdRmewMIWpZyZWhLgEt921LzT9I0peySF-5X_44EqSgcNZWpVrf2DKUD_6MSXvu1d6Rjgl-ohXzTRru9ZK6N7gNx6iVczXY7PZ6j5y-d0Y_6CURHx8q8uM6ETlLt-TuGmJo15kosve81XvrpjxHz8oaxqayEzqM2xmCg87eQJz_W5UC2b7ms6IOH_ZuRnFggtv_dmTe2JXwaFCepu7IlTuHp1MaI_RKsgE7zmEQ7d2TjRm2J4dCJc7lHSBYAKBgGruo6pxZIaxcCRqaAATvdCV26I8JNYySjDPwL4Gwx6Zg-zniUyyEg1iUhBMA2sbwtu8LHnJaVw5ERWC9IIPzVq2SjO7ceIJ8DTvuUnLjJ5PzqXQUqEXK0Ulrfpa16oot4WO8ueAFvZbZXm1N05KoayGOpIP4xz1V0DtjfLLpUAeKEz_K-gjw5bV2hbhuSSs7ZByOo87j8XX4hlpWFGr1TD0jYVTZgn1E5k-EAdraHthH3u-3A--IPtX3Jp05IxLrrIbp-oiJBCNiayLqSoVleyP3xjk9o4BK9frge6XFpnKaNZ48RklVQZcvZjEcPtkEWiUKii-tqxblrS8Kz8zPjvARoAW1NmeYEzJg1SXBMzwcS3257rfGWCItT8twrO2XLCPyvamuQJ9UtmrOhaRtg7hpU9E1dSzYy3spCnsgSxgZAr726T1xiWuw4wGZKdu01zRlByj66QjIpNLAL4er4qJ5lH9lIPlhAVoJAAVK6aoPD-J09I3kNYk58TgR92SxhHiRF5Ihg4jsZMpbwpEjB5r4ean1sKBDS7sQkrNuMchfQyz9gQfSNahp5AfK5OghHfX_6X0L9aNTRrFXcvbRrhEj13_iHekGy5l4BfzRcFkerrtaxRla1TAlFyIRAemfSQ40XjKZiUueArDjUGeAuDyJc2zazFYUgot_Rm7aAR_ihQ.yX7tIFnrHKcoeXqL7r3EVA"
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