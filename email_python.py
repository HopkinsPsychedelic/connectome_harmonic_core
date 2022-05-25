#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 15 13:02:23 2022

@author: bwinston
"""
import yagmail 

user = 'chapnotifications@gmail.com'
app_password = 'bgzflqzszjjdusip' # a token for gmail

def send_email_notif(recipient = 'chapnotifications@gmail.com',subject = 'notification',content = 'chap finished'):
    with yagmail.SMTP(user, app_password) as yag:
        yag.send(recipient, subject, content)